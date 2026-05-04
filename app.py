import gzip
import json
import textwrap
import webbrowser
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter import LEFT, W, E, S, N, RAISED, Toplevel, Label, Button

import pandas as pd
import requests
import graphviz


# ── Constants ────────────────────────────────────────────────────────────────

# Accepts either the compressed (.json.gz) or plain (.json) data file
_BASE = Path(__file__).parent
DATA_FILE = next(
    (p for p in [
        _BASE / "ORTHOLOGY-ALLIANCE-JSON_COMBINED.json.gz",
        _BASE / "ORTHOLOGY-ALLIANCE-JSON_COMBINED.json",
        _BASE / "ORTHOLOGY-ALLIANCE-JSON_COMBINED_37.json.gz",
        _BASE / "ORTHOLOGY-ALLIANCE-JSON_COMBINED_37.json",
    ] if p.exists()),
    None,
)
API_BASE   = "https://www.alliancegenome.org/api/gene/"
GENE_PAGE  = "https://www.alliancegenome.org/gene/"
AMIGO_URL  = "https://amigo.geneontology.org/amigo/term/"
MYGENE_URL = "https://mygene.info/v3/query"

# MyGene.info species names for the input (zebrafish) genes
MYGENE_SPECIES = {
    "Danio rerio":             "zebrafish",
    "Homo sapiens":            "human",
    "Mus musculus":            "mouse",
    "Rattus norvegicus":       "rat",
    "Drosophila":              "fruitfly",
    "Caenorhabditis elegans":  "nematode",
    "Saccharomyces cerevisiae":"yeast",
}

ASPECT_KEYS    = ["BP", "MF", "CC"]
ASPECT_LABELS  = {"BP": "Biological Process", "MF": "Molecular Function", "CC": "Cellular Component"}
ASPECT_COLORS  = {"BP": "palegreen",          "MF": "lightyellow",        "CC": "lightcoral"}

OUTPUT_SPECIES = [
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Drosophila",
    "Caenorhabditis elegans",
    "Saccharomyces cerevisiae",
]
COUNT_OPTIONS = [str(i) for i in range(1, 12)]


# ── Helpers ──────────────────────────────────────────────────────────────────

def wrap_text(text: str, width: int = 100) -> str:
    return "\n".join(textwrap.wrap(str(text), width))


def fetch_go_annotations(gene_symbol: str, species: str = "zebrafish") -> dict:
    """Return GO data from MyGene.info keyed by aspect (BP, MF, CC).
    Each value is a list of dicts with keys: id, term, evidence, qualifier."""
    try:
        r = requests.get(
            MYGENE_URL,
            params={"q": gene_symbol, "species": species, "fields": "go", "size": 1},
            timeout=10,
        )
        r.raise_for_status()
        hits = r.json().get("hits", [])
        if not hits:
            return {}
        go = hits[0].get("go", {})
        # Normalise: each aspect value may be a single dict instead of a list
        return {k: (v if isinstance(v, list) else [v]) for k, v in go.items() if k in ASPECT_KEYS}
    except requests.RequestException:
        return {}


def build_go_graph(gene_go_data: dict[str, dict], title: str = "go_annotations") -> graphviz.Digraph:
    """Build a Graphviz diagram for one or more genes' GO annotations.

    Layout: gene node → aspect proxy node → GO term nodes.
    Aspect nodes are colour-coded; one subgraph cluster per gene.
    """
    g = graphviz.Digraph(title, filename=f"{title}.gv")
    g.attr(rankdir="LR", fontsize="11")

    for gene_symbol, go_data in gene_go_data.items():
        with g.subgraph(name=f"cluster_{gene_symbol}") as cluster:
            cluster.attr(label=gene_symbol, style="rounded,filled",
                         fillcolor="lightblue2", color="steelblue")
            cluster.node(
                gene_symbol, label=gene_symbol,
                shape="box", style="filled", fillcolor="lightblue2", fontsize="12",
            )
            for aspect in ASPECT_KEYS:
                terms = go_data.get(aspect, [])
                if not terms:
                    continue
                aspect_node = f"{gene_symbol}_{aspect}"
                cluster.node(
                    aspect_node, label=ASPECT_LABELS[aspect],
                    shape="box", style="filled",
                    fillcolor=ASPECT_COLORS[aspect], fontsize="10",
                )
                g.edge(gene_symbol, aspect_node)
                for term in terms:
                    go_id   = term.get("id", "")
                    go_name = term.get("term", "")
                    term_node = f"{gene_symbol}_{go_id}"
                    label = wrap_text(f"{go_id}\n{go_name}", 30)
                    cluster.node(
                        term_node, label=label,
                        shape="box", style="filled",
                        fillcolor=ASPECT_COLORS[aspect], fontsize="9",
                    )
                    g.edge(aspect_node, term_node)
    return g


def fetch_description(gene_id: str) -> str:
    """Fetch automated gene synopsis from Alliance Genome API."""
    try:
        response = requests.get(API_BASE + gene_id, timeout=10)
        response.raise_for_status()
        return response.json().get("automatedGeneSynopsis", "Description not found")
    except requests.RequestException:
        return "Description not found"


# ── GO Annotation Window ──────────────────────────────────────────────────────

class GOAnnotationWindow:
    """Popup that shows GO annotations for one or more genes.

    gene_go_data: {gene_symbol: {"BP": [...], "MF": [...], "CC": [...]}}
    Each term dict has keys: id, term, evidence, qualifier.
    """

    def __init__(self, parent: tk.Tk, gene_go_data: dict[str, dict]):
        self.parent = parent
        self.gene_go_data = gene_go_data
        self.genes = list(gene_go_data.keys())

        self.win = Toplevel(parent)
        self.win.title("GO Annotations")
        self.win.geometry("900x600")
        self.win.transient(parent)

        self._build_ui()
        if self.genes:
            self._load_gene(self.genes[0])

    def _build_ui(self):
        win = self.win

        # ── Top bar: gene selector ────────────────────────────────────────────
        top = tk.Frame(win)
        top.pack(fill=tk.X, padx=8, pady=6)

        Label(top, text="Gene:", font=("Helvetica", 10, "bold")).pack(side=tk.LEFT)
        self.gene_var = tk.StringVar()
        selector = ttk.Combobox(
            top, textvariable=self.gene_var,
            values=self.genes, state="readonly", width=20,
        )
        selector.pack(side=tk.LEFT, padx=6)
        selector.bind("<<ComboboxSelected>>", lambda _: self._load_gene(self.gene_var.get()))

        self.count_label = Label(top, text="", fg="gray")
        self.count_label.pack(side=tk.LEFT, padx=10)

        # ── Notebook: one tab per aspect ──────────────────────────────────────
        self.notebook = ttk.Notebook(win)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        self.trees: dict[str, ttk.Treeview] = {}
        for aspect in ASPECT_KEYS:
            frame = tk.Frame(self.notebook)
            self.notebook.add(frame, text=ASPECT_LABELS[aspect])

            cols = ("go_id", "term", "evidence", "qualifier")
            tree = ttk.Treeview(frame, columns=cols, show="headings", height=18)
            tree.heading("go_id",     text="GO ID",      command=lambda c="go_id",     t=aspect: self._sort(t, c))
            tree.heading("term",      text="Term",        command=lambda c="term",      t=aspect: self._sort(t, c))
            tree.heading("evidence",  text="Evidence",    command=lambda c="evidence",  t=aspect: self._sort(t, c))
            tree.heading("qualifier", text="Qualifier",   command=lambda c="qualifier", t=aspect: self._sort(t, c))
            tree.column("go_id",     width=110, anchor=W)
            tree.column("term",      width=420, anchor=W)
            tree.column("evidence",  width=90,  anchor=W)
            tree.column("qualifier", width=180, anchor=W)
            tree.tag_configure("go_row", background=ASPECT_COLORS[aspect])

            sb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
            tree.configure(yscrollcommand=sb.set)
            tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            sb.pack(side=tk.RIGHT, fill=tk.Y)

            tree.bind("<Double-1>", self._on_term_click)
            self.trees[aspect] = tree

        # ── Bottom buttons ────────────────────────────────────────────────────
        btn_frame = tk.Frame(win)
        btn_frame.pack(fill=tk.X, padx=8, pady=6)

        Button(btn_frame, text="Graph this gene",
               command=self._graph_current).pack(side=tk.LEFT, padx=4)
        Button(btn_frame, text="Graph all genes",
               command=self._graph_all).pack(side=tk.LEFT, padx=4)
        Button(btn_frame, text="Close",
               command=self.win.destroy).pack(side=tk.RIGHT, padx=4)

    # ── Loading ───────────────────────────────────────────────────────────────

    def _load_gene(self, gene_symbol: str):
        self.gene_var.set(gene_symbol)
        go_data = self.gene_go_data.get(gene_symbol, {})
        total = sum(len(go_data.get(a, [])) for a in ASPECT_KEYS)
        self.count_label.config(text=f"{total} annotation(s)")

        for aspect in ASPECT_KEYS:
            tree = self.trees[aspect]
            tree.delete(*tree.get_children())
            for term in go_data.get(aspect, []):
                tree.insert("", tk.END, values=(
                    term.get("id", ""),
                    term.get("term", ""),
                    term.get("evidence", ""),
                    term.get("qualifier", ""),
                ), tags=("go_row",))

    # ── Sorting ───────────────────────────────────────────────────────────────

    def _sort(self, aspect: str, col: str):
        tree = self.trees[aspect]
        rows = [(tree.set(iid, col), iid) for iid in tree.get_children()]
        rows.sort(key=lambda x: x[0].lower())
        for i, (_, iid) in enumerate(rows):
            tree.move(iid, "", i)

    # ── Double-click GO term → AmiGO2 ────────────────────────────────────────

    def _on_term_click(self, event):
        widget = event.widget
        selection = widget.selection()
        if not selection:
            return
        go_id = widget.set(selection[0], "go_id")  # e.g. "GO:0001234"
        if go_id.startswith("GO:"):
            webbrowser.open(AMIGO_URL + go_id)

    # ── Graphing ──────────────────────────────────────────────────────────────

    def _graph_current(self):
        gene = self.gene_var.get()
        if not gene:
            return
        g = build_go_graph({gene: self.gene_go_data[gene]}, title=f"go_{gene}")
        g.view()

    def _graph_all(self):
        g = build_go_graph(self.gene_go_data, title="go_all_genes")
        g.view()


# ── Main Application ─────────────────────────────────────────────────────────

class GeneSearchApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Gene Orthology Search")
        self.root.attributes("-fullscreen", True)

        # State
        self.desc_cache: dict[str, str] = {}
        self.graph: graphviz.Digraph | None = None
        self._last_searched_genes: list[str] = []   # gene symbols from most recent search

        # Filter variables
        self.species_vars: dict[str, tk.BooleanVar] = {
            s: tk.BooleanVar() for s in OUTPUT_SPECIES
        }
        self.count_vars: dict[str, tk.BooleanVar] = {
            c: tk.BooleanVar() for c in COUNT_OPTIONS
        }
        self.var_all_species = tk.BooleanVar()
        self.var_all_count = tk.BooleanVar()

        # Column visibility
        self.var_show_orthology = tk.IntVar(value=1)
        self.var_show_description = tk.IntVar(value=1)

        self._build_ui()
        self._load_data()

    # ── Data loading ──────────────────────────────────────────────────────────

    def _load_data(self):
        if DATA_FILE is None:
            messagebox.showerror(
                "Data file missing",
                "Could not find ORTHOLOGY-ALLIANCE-JSON_COMBINED.json(.gz) next to app.py.\n"
                "Download it from https://www.alliancegenome.org/downloads and place it here:\n"
                f"{_BASE}",
            )
            self._raw_data = []
            return
        try:
            opener = gzip.open if DATA_FILE.suffix == ".gz" else open
            with opener(DATA_FILE) as f:
                self._raw_data = json.load(f)["data"]
        except (json.JSONDecodeError, KeyError) as e:
            messagebox.showerror("Data file error", f"Failed to parse JSON:\n{e}")
            self._raw_data = []

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        root = self.root
        root.columnconfigure(1, weight=1)
        root.rowconfigure(4, weight=1)

        # Row 0: gene ID entry
        Label(root, text="Enter Gene IDs (semicolon-separated)", font=("Helvetica", 10)).grid(
            row=0, column=0, sticky=W, pady=4, padx=5
        )
        self.gene_entry = tk.Entry(root)
        self.gene_entry.grid(row=0, column=1, sticky=W + E, pady=4, padx=5)
        self._build_species_menu(row=0, col=2)

        # Row 1: file upload
        Button(root, text="Choose Excel file", command=self._upload_file).grid(
            row=1, column=0, sticky=W + E, pady=4, padx=5
        )
        self.file_entry = tk.Entry(root)
        self.file_entry.grid(row=1, column=1, sticky=W + E, pady=4, padx=5)
        self._build_count_menu(row=1, col=2)

        # Row 2: search button
        Button(root, text="Search", command=self._search).grid(
            row=2, column=2, sticky=W + E, pady=4, padx=5
        )

        # Row 3: column toggles + GO button
        Label(root, text="Display columns", font=("Helvetica", 10)).grid(
            row=3, column=0, sticky=W, pady=4, padx=5
        )
        tk.Checkbutton(
            root, text="Orthology", variable=self.var_show_orthology,
            command=self._update_columns,
        ).grid(row=3, column=0, sticky=W, padx=120)
        tk.Checkbutton(
            root, text="Description", variable=self.var_show_description,
            command=self._update_columns,
        ).grid(row=3, column=1, sticky=W)
        self.go_btn = Button(
            root, text="GO Annotation",
            state=tk.DISABLED, command=self._on_go_annotation_click,
        )
        self.go_btn.grid(row=3, column=2, sticky=W + E, pady=4, padx=5)

        # Row 4: results tree
        self._build_tree(row=4)

        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        Label(root, textvariable=self.status_var, anchor=W, relief=tk.SUNKEN).grid(
            row=5, column=0, columnspan=3, sticky=W + E, padx=5, pady=2
        )

    def _build_species_menu(self, row: int, col: int):
        mb = tk.Menubutton(self.root, text="Output species", relief=RAISED)
        mb.menu = tk.Menu(mb, tearoff=0)
        mb["menu"] = mb.menu
        mb.menu.add_checkbutton(
            label="Select all", variable=self.var_all_species,
            command=self._toggle_all_species,
        )
        mb.menu.add_separator()
        for species in OUTPUT_SPECIES:
            mb.menu.add_checkbutton(
                label=species, variable=self.species_vars[species],
                command=lambda: self.var_all_species.set(False),
            )
        mb.grid(row=row, column=col, sticky=E + W, pady=4, padx=5)

    def _build_count_menu(self, row: int, col: int):
        mb = tk.Menubutton(self.root, text="Algorithm count", relief=RAISED)
        mb.menu = tk.Menu(mb, tearoff=0)
        mb["menu"] = mb.menu
        mb.menu.add_checkbutton(
            label="Select all", variable=self.var_all_count,
            command=self._toggle_all_count,
        )
        mb.menu.add_separator()
        for count in COUNT_OPTIONS:
            mb.menu.add_checkbutton(
                label=count, variable=self.count_vars[count],
                command=lambda: self.var_all_count.set(False),
            )
        mb.grid(row=row, column=col, sticky=E + W, pady=4, padx=5)

    def _build_tree(self, row: int):
        ttk.Style().configure("Treeview", rowheight=45)
        cols = ("outputSpecies", "outputGeneSymbol", "Count", "Methods", "Description", "Link")
        self.tree = ttk.Treeview(self.root, columns=cols)

        self.tree.column("#0", width=120)
        self.tree.heading("#0", text="Input Gene")
        self.tree.column("outputSpecies", anchor=W, width=160)
        self.tree.heading("outputSpecies", text="Output Species")
        self.tree.column("outputGeneSymbol", anchor=W, width=100)
        self.tree.heading("outputGeneSymbol", text="Output Gene")
        self.tree.column("Count", anchor=W, width=50)
        self.tree.heading("Count", text="Count")
        self.tree.column("Methods", anchor=W, width=200)
        self.tree.heading("Methods", text="Orthology")
        self.tree.column("Description", anchor=W, width=450)
        self.tree.heading("Description", text="Description")
        self.tree.column("Link", width=100)
        self.tree.heading("Link", text="Reference")

        self.tree.bind("<Double-1>", self._on_double_click)
        self.tree.grid(row=row, column=0, columnspan=3, sticky=S + N + W + E, pady=4, padx=5)

        scrollbar = ttk.Scrollbar(self.root, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        scrollbar.grid(row=row, column=3, sticky=S + N)

    # ── Filter toggle helpers ─────────────────────────────────────────────────

    def _toggle_all_species(self):
        val = self.var_all_species.get()
        for var in self.species_vars.values():
            var.set(val)

    def _toggle_all_count(self):
        val = self.var_all_count.get()
        for var in self.count_vars.values():
            var.set(val)

    def _selected_species(self) -> list[str]:
        return [s for s, v in self.species_vars.items() if v.get()]

    def _selected_counts(self) -> set[str]:
        return {c for c, v in self.count_vars.items() if v.get()}

    # ── Column visibility ─────────────────────────────────────────────────────

    def _update_columns(self):
        show_orth = bool(self.var_show_orthology.get())
        show_desc = bool(self.var_show_description.get())
        cols = ["outputSpecies", "outputGeneSymbol", "Count"]
        if show_orth:
            cols.append("Methods")
        if show_desc:
            cols.append("Description")
        cols.append("Link")
        self.tree["displaycolumns"] = cols

    # ── Gene name extraction ──────────────────────────────────────────────────

    def _get_gene_names(self) -> list[str]:
        """Return gene names from the text entry, or from the uploaded Excel file."""
        raw = self.gene_entry.get().strip()
        if raw:
            return [n.strip() for n in raw.split(";") if n.strip()]

        filename = self.file_entry.get().strip()
        if not filename:
            messagebox.showwarning("No input", "Enter gene IDs or choose an Excel file.")
            return []
        try:
            df = pd.read_excel(filename)
            return [str(v).strip() for v in df.iloc[:, 1] if str(v).strip()]
        except Exception as e:
            messagebox.showerror("File error", f"Could not read Excel file:\n{e}")
            return []

    # ── Search ────────────────────────────────────────────────────────────────

    def _search(self):
        gene_names = self._get_gene_names()
        if not gene_names:
            return

        selected_species = self._selected_species()
        selected_counts = self._selected_counts()

        if not selected_species:
            messagebox.showwarning("No species", "Select at least one output species.")
            return
        if not selected_counts:
            messagebox.showwarning("No count", "Select at least one algorithm count.")
            return

        # Filter the raw data
        gene_name_set = set(gene_names)
        matches = [
            s for s in self._raw_data
            if s["Gene2SpeciesName"] in selected_species
            and s["Gene1Symbol"] in gene_name_set
            and str(s["AlgorithmsMatch"]) in selected_counts
        ]

        if not matches:
            messagebox.showinfo("No results", "No ortholog matches found for the given filters.")
            return

        self.tree.delete(*self.tree.get_children())
        self.desc_cache.clear()
        self.graph = self._new_graph()

        self.status_var.set(f"Fetching descriptions for {len(matches)} results…")
        self.root.update_idletasks()

        # Fetch all descriptions concurrently
        unique_gene1_ids = list({s["Gene1ID"] for s in matches})
        descriptions: dict[str, str] = {}
        with ThreadPoolExecutor(max_workers=10) as pool:
            future_to_id = {pool.submit(fetch_description, gid): gid for gid in unique_gene1_ids}
            for future in as_completed(future_to_id):
                gid = future_to_id[future]
                descriptions[gid] = future.result()

        # Populate tree
        seen_genes: set[str] = set()
        for entry in matches:
            gene_symbol = entry["Gene1Symbol"]
            gene_id = entry["Gene1ID"]
            out_species = entry["Gene2SpeciesName"]
            out_symbol = entry["Gene2Symbol"]
            count = str(entry["AlgorithmsMatch"])
            methods = entry.get("Algorithms", [])
            full_desc = descriptions.get(gene_id, "Description not found")
            link = GENE_PAGE + gene_id

            self._add_graph_edges(gene_symbol, out_species, out_symbol, methods)

            first_occurrence = gene_symbol not in seen_genes
            if first_occurrence:
                seen_genes.add(gene_symbol)
                short_desc = full_desc[:60] + "… (double-click for full)" if len(full_desc) > 60 else full_desc
                self.desc_cache[gene_symbol] = full_desc
                self.tree.insert("", tk.END, iid=gene_symbol, text=gene_symbol)
                display_desc = wrap_text(short_desc, 60)
                display_link = link
            else:
                display_desc = ""
                display_link = ""

            self.tree.insert(
                gene_symbol, tk.END,
                values=(
                    out_species,
                    out_symbol,
                    count,
                    wrap_text("-".join(methods), 45),
                    display_desc,
                    display_link,
                ),
            )

        self.graph.view()
        self._last_searched_genes = list(seen_genes)
        self.go_btn.config(state=tk.NORMAL)
        self.status_var.set(f"Done — {len(matches)} results for {len(seen_genes)} gene(s).")

    # ── Graph ─────────────────────────────────────────────────────────────────

    def _new_graph(self) -> graphviz.Digraph:
        g = graphviz.Digraph(
            "orthology",
            filename="orthology.gv",
            strict=True,
            node_attr={"color": "lightblue2", "style": "filled"},
        )
        g.attr("node", shape="box")
        g.node("gene", label="Gene ID")
        return g

    def _add_graph_edges(self, gene_symbol: str, out_species: str, out_symbol: str, methods: list):
        assert self.graph is not None
        self.graph.node(gene_symbol)
        self.graph.edge("gene", gene_symbol)
        self.graph.edge(gene_symbol, out_species)
        self.graph.edge(out_species, out_symbol)
        for method in methods:
            self.graph.edge(out_symbol, method)

    # ── Tree interactions ─────────────────────────────────────────────────────

    def _on_double_click(self, event):
        selection = self.tree.selection()
        if not selection:
            return
        item = self.tree.item(selection[0])
        values = item.get("values", [])

        # Top-level gene row (iid == gene symbol) → open description window
        gene_symbol = selection[0]
        if gene_symbol in self.desc_cache:
            self._open_description_window(gene_symbol)
        else:
            # Child row — open link
            link = values[-1] if values else ""
            if link:
                webbrowser.open(link)

    def _open_description_window(self, gene_symbol: str):
        full_desc = self.desc_cache.get(gene_symbol, "")
        if not full_desc:
            return
        win = Toplevel(self.root)
        win.title(f"Description — {gene_symbol}")
        win.transient(self.root)
        Label(win, text=wrap_text(full_desc, 80), justify=LEFT, wraplength=600, padx=10, pady=10).pack()
        link = GENE_PAGE + gene_symbol
        Button(win, text="Open reference page", command=lambda: webbrowser.open(link)).pack(pady=4)
        Button(win, text="Close", command=win.destroy).pack(pady=4)

    # ── GO Annotation ────────────────────────────────────────────────────────

    def _on_go_annotation_click(self):
        if not self._last_searched_genes:
            messagebox.showinfo("No genes", "Run a search first.")
            return

        self.status_var.set(f"Fetching GO annotations for {len(self._last_searched_genes)} gene(s)…")
        self.root.update_idletasks()

        gene_go_data: dict[str, dict] = {}
        with ThreadPoolExecutor(max_workers=10) as pool:
            futures = {
                pool.submit(fetch_go_annotations, sym): sym
                for sym in self._last_searched_genes
            }
            for future in as_completed(futures):
                sym = futures[future]
                gene_go_data[sym] = future.result()

        empty = [g for g, d in gene_go_data.items() if not d]
        if empty:
            self.status_var.set(
                f"GO annotations loaded (no data for: {', '.join(empty)})."
            )
        else:
            self.status_var.set("GO annotations loaded.")

        GOAnnotationWindow(self.root, gene_go_data)

    # ── File upload ───────────────────────────────────────────────────────────

    def _upload_file(self):
        filename = filedialog.askopenfilename(
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if filename:
            self.file_entry.delete(0, tk.END)
            self.file_entry.insert(0, filename)


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    root = tk.Tk()
    app = GeneSearchApp(root)
    root.mainloop()
