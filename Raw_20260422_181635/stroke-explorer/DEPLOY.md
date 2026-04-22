# Stroke Cell Death Explorer — GitHub Pages Deployment

## Folder structure
```
stroke-explorer/
├── index.html          (main app, ~44 KB)
└── data/
    ├── datasets.json   (4.8 KB)  — 9 dataset metadata cards
    ├── heatmap.json    (28 KB)   — cell death z-scores per cell type × dataset
    ├── forest.json     (6.5 KB)  — weighted meta-analysis results
    ├── temporal.json   (32 KB)   — time course data (GSE197341 + GSE279666)
    ├── glial.json      (42 KB)   — reactive glia state scores
    └── gene_expr.json  (7.2 MB)  — per-gene expression across all 9 datasets
```

## Deploy to GitHub Pages (from scratch)

### Step 1: Create a new GitHub repository
Go to https://github.com/new and create a new public repository (e.g. `stroke-cell-death-explorer`).

### Step 2: Initialize and push

```bash
# Navigate to the stroke-explorer folder
cd /path/to/stroke-explorer

# Initialize git
git init
git add .
git commit -m "Initial commit: Stroke Cell Death Explorer"

# Add your GitHub remote (replace USERNAME and REPO)
git remote add origin https://github.com/USERNAME/stroke-cell-death-explorer.git

# Push
git branch -M main
git push -u origin main
```

### Step 3: Enable GitHub Pages
1. Go to your repository on GitHub
2. Click **Settings** → **Pages** (left sidebar)
3. Under **Source**, select **Deploy from a branch**
4. Branch: `main`, Folder: `/ (root)`
5. Click **Save**

Your app will be live at:
`https://USERNAME.github.io/stroke-cell-death-explorer/`

(GitHub Pages deployment takes ~1-2 minutes after the first push.)

## Notes
- The `gene_expr.json` file is 7.2 MB. GitHub Pages serves it fine as a static file.
- The app loads all data on first visit; subsequent tab switches are instant.
- No server-side code required — fully static.
- Works offline once loaded in the browser.
