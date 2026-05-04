# Gene Orthology & GO Annotation Explorer

A desktop application for exploring orthologous genes across model organisms and visualizing Gene Ontology (GO) annotations. Built on top of the [Alliance of Genome Resources](https://www.alliancegenome.org/) orthology dataset and the [MyGene.info](https://mygene.info/) annotation API.

---

## Features

- **Orthology search** — query by gene symbol (or upload an Excel file) and filter results by output species and number of supporting algorithms
- **Concurrent API calls** — gene descriptions and GO annotations are fetched in parallel, keeping the interface responsive
- **Interactive results table** — expandable rows show orthologs per gene; toggle Orthology and Description columns on/off
- **GO Annotation viewer** — tabbed window (Biological Process · Molecular Function · Cellular Component) with sortable columns; double-click any GO term to open it in AmiGO2
- **Graphviz diagrams** — one-click generation of orthology relationship graphs and GO annotation graphs (single gene or all genes at once)

---

## Screenshots

> _Add screenshots here once the app is running._

---

## Requirements

| Requirement | Version |
|---|---|
| Python | 3.10 or higher |
| Graphviz (system) | any recent version |

Install the Graphviz system tool before running the app:

```bash
# macOS
brew install graphviz

# Ubuntu / Debian
sudo apt install graphviz

# Windows
# Download the installer from https://graphviz.org/download/
```

---

## Setup

```bash
# 1. Clone the repository
git clone <your-repo-url>
cd gene-orthology-explorer

# 2. Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate       # macOS / Linux
.venv\Scripts\activate          # Windows

# 3. Install Python dependencies
pip install -r requirements.txt

# 4. Download the orthology data file
#    Go to: https://www.alliancegenome.org/downloads  →  Orthology section
#    Download: ORTHOLOGY-ALLIANCE-JSON_COMBINED.json.gz
#    Place it in the same folder as app.py

# 5. Run
python app.py
```

The app accepts both the compressed (`.json.gz`) and uncompressed (`.json`) versions of the data file — no manual decompression needed.

---

## Usage

### Orthology Search

1. **Enter gene IDs** in the text field, separated by semicolons — e.g. `sox2;pax6;tp53`  
   **or** click **Choose Excel file** to load a spreadsheet (gene symbols must be in column B).
2. Click **Output species** and select one or more target species.
3. Click **Algorithm count** and select which confidence levels to include (higher = more algorithms agree).
4. Click **Search**. The status bar tracks progress while descriptions are fetched from the Alliance Genome API.
5. Results appear grouped by input gene. Click the arrow next to a gene to expand its ortholog rows.
6. **Double-click** a gene row to read its full description and open its Alliance Genome reference page.

### GO Annotations

1. Run a search first — the **GO Annotation** button activates once results are loaded.
2. Click **GO Annotation** to open the annotation viewer. Annotations are fetched for all genes in the search simultaneously.
3. Use the gene selector dropdown to switch between genes.
4. Click any column header to sort that tab's results.
5. **Double-click** a GO term row to open that term in [AmiGO2](https://amigo.geneontology.org/).

### Graphs

| Button | What it generates |
|---|---|
| _(auto, after Search)_ | Orthology diagram: input gene → output species → ortholog → supporting methods |
| **Graph this gene** | GO annotation diagram for the currently selected gene |
| **Graph all genes** | GO annotation diagram for every gene in the search, in one view |

Graph files are saved as `.gv` and rendered as PDF in the same folder as `app.py`.  
GO term nodes are colour-coded: **green** = Biological Process · **yellow** = Molecular Function · **red** = Cellular Component.

---

## Data Sources

| Source | Used for |
|---|---|
| [Alliance of Genome Resources](https://www.alliancegenome.org/) | Orthology relationships and gene descriptions |
| [MyGene.info](https://mygene.info/) | GO annotations (BP, MF, CC) |
| [AmiGO2](https://amigo.geneontology.org/) | GO term reference pages (opened on double-click) |

---

## Project Structure

```
.
├── app.py             # Application — search, GO annotation, graphing
├── requirements.txt   # Python dependencies
└── README.md
```

---

## Supported Species

**Output species** (ortholog targets):

- *Homo sapiens*
- *Mus musculus*
- *Rattus norvegicus*
- *Drosophila melanogaster*
- *Caenorhabditis elegans*
- *Saccharomyces cerevisiae*

Input genes are assumed to be *Danio rerio* (zebrafish) symbols, consistent with the Alliance orthology dataset.

---

## Acknowledgements

Developed as part of CS 505 – Special Topics (Spring 2021) under Dr. Jamil at the Computer Science Department at the University of Idaho. 
Orthology data provided by the [Alliance of Genome Resources](https://www.alliancegenome.org/).  
GO annotation data provided by [MyGene.info](https://mygene.info/).
