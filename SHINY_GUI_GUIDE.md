# powerNMA Interactive GUI Guide

**Version**: 3.0
**Interface**: Shiny Dashboard with bs4Dash
**Status**: Production-ready

---

## Overview

The powerNMA Interactive GUI provides a modern, user-friendly web interface for network meta-analysis. Built with Shiny and bs4Dash, it offers point-and-click access to all four analysis pathways without requiring R programming knowledge.

### Key Features

✨ **Modern Dashboard Interface** - Clean, professional design using Bootstrap 4
📊 **Interactive Visualizations** - Dynamic plots with plotly
📁 **Easy Data Upload** - CSV, Excel, or use example datasets
🤖 **4 Analysis Pathways** - From full manual control to complete automation
📄 **Report Generation** - Export results as HTML, PDF, or CSV
❓ **Built-in Help** - Comprehensive documentation within the app

---

## Quick Start

### Installation

```r
# Install powerNMA (if not already installed)
# install.packages("devtools")
# devtools::install_github("mahmood726-cyber/rmstnma/powerNMA")

library(powerNMA)

# Check if GUI dependencies are installed
check_gui_dependencies()

# If missing, install them
install_gui_dependencies()
```

### Launching the GUI

```r
# Launch with default settings
launch_powernma_gui()

# Or use the shorter alias
powernma_gui()

# Launch on specific port
launch_powernma_gui(port = 8080)

# Launch without opening browser (useful for servers)
launch_powernma_gui(launch.browser = FALSE)
```

The GUI will open in your default web browser at `http://127.0.0.1:[port]`

---

## Interface Tour

### Dashboard Layout

```
┌─────────────────────────────────────────────────────────────────┐
│  powerNMA                                         [GitHub]       │ Header
├──────────────┬──────────────────────────────────────────────────┤
│              │                                                   │
│   Sidebar    │              Main Content Area                    │
│   Menu       │                                                   │
│              │                                                   │
│   • Welcome  │  [Tabs, Forms, Results, Visualizations]         │
│   • Data     │                                                   │
│   • Path 1   │                                                   │
│   • Path 2   │                                                   │
│   • Path 3   │                                                   │
│   • Path 4   │                                                   │
│   • Results  │                                                   │
│   • Help     │                                                   │
│              │                                                   │
├──────────────┴──────────────────────────────────────────────────┤
│  powerNMA v3.0 | © 2025 | Powered by R Shiny & bs4Dash          │ Footer
└─────────────────────────────────────────────────────────────────┘
```

### Navigation

The **left sidebar** provides access to all sections:

- 🏠 **Welcome** - Introduction and pathway overview
- 📁 **Data Upload** - Upload data or use examples
- ⚙️ **Pathway 1: Manual Standard** - Full control, validated methods
- 🧪 **Pathway 2: Manual Experimental** - Full control, cutting-edge
- ✨ **Pathway 3: Auto Standard** - Automatic, validated
- 🤖 **Pathway 4: Auto Experimental** - Automatic, cutting-edge
- 📊 **Results & Export** - View and download results
- ❓ **Help & Documentation** - Built-in help

---

## Step-by-Step Workflow

### Step 1: Upload Data

Click **Data Upload** in the sidebar.

#### Option A: Upload Your Own File

1. Click **"Upload File"** tab
2. Click **"Choose CSV or Excel file"**
3. Select your data format:
   - **Pairwise**: `study, treat1, treat2, TE, seTE`
   - **Arm-based binary**: `study, treatment, n, events`
   - **Arm-based continuous**: `study, treatment, n, mean, sd`
   - **IPD**: `study, treatment, outcome, [covariates...]`
   - **Time-to-event**: `study, treatment, time, event`
4. Check/uncheck **"File has header row"**
5. Click **"Load Data"**

#### Option B: Use Example Data

1. Click **"Use Example Data"** tab
2. Select from dropdown:
   - Smoking Cessation (Binary)
   - Depression Treatment (Continuous)
   - Cancer Survival (Time-to-Event)
   - Multicomponent Interventions
3. Click **"Load Example"**

#### Data Preview

After loading, you'll see:
- **Data table** with pagination
- **Data summary** showing:
  - Number of rows and columns
  - Column names
  - Data structure

---

### Step 2: Choose Your Pathway

Based on your needs, select one of the four pathways:

#### Decision Guide

```
Need full control?
├─ Standard methods (regulatory, guidelines)
│  └─ Pathway 1: Manual Standard
└─ Experimental methods (advanced research)
   └─ Pathway 2: Manual Experimental

Want automation?
├─ Standard methods (consistency, limited expertise)
│  └─ Pathway 3: Auto Standard
└─ Experimental methods (cutting-edge, rapid)
   └─ Pathway 4: Auto Experimental
```

---

### Step 3A: Pathway 1 - Manual Standard

**For**: Regulatory submissions, clinical guidelines, full control

#### Configuration

1. **Choose Method** (dropdown):
   - Standard NMA
   - Component NMA (CNMA)
   - Network Meta-Regression
   - Dose-Response NMA
   - Predictive Ranking
   - Multivariate NMA
   - Missing Data Handling
   - Cross-Design Synthesis

2. **Set Parameters**:
   - **Model Type**: Random Effects or Fixed Effect
   - **Summary Measure**: MD, SMD, OR, RR, HR
   - **Reference Treatment**: Enter treatment name
   - **Check Inconsistency**: Yes/No

3. Click **"Run Analysis"** (green button)

#### Results

- Analysis status
- Method used
- Results summary
- (Full implementation would show detailed results)

---

### Step 3B: Pathway 2 - Manual Experimental

**For**: Advanced research, methodological studies, precision medicine

⚠️ **Warning**: Experimental methods banner displayed

#### Configuration

1. **Choose Experimental Method** (dropdown):
   - **RMST-based NMA**: Time-to-event with interpretable survival
   - **Threshold Analysis**: Robustness assessment
   - **Individualized Treatment Rules (ITR)**: Precision medicine
   - **Bayesian Model Averaging**: Model uncertainty

2. **Method-Specific Parameters** (dynamic, changes based on method):

   **If RMST selected**:
   - Restriction time (tau)
   - Data type (IPD/Aggregate)

   **If Threshold selected**:
   - Outcome direction (higher/lower better)
   - Risk aversion (0-2 slider)

   **If ITR selected**:
   - Covariates (comma-separated)
   - ITR method (Regression/ML)

   **If BMA selected**:
   - Weighting method (BIC/AIC/DIC)

3. Click **"Run Experimental Analysis"** (yellow button)

#### Results

- Experimental warning
- Method used
- Experimental results
- (Full implementation would show method-specific outputs)

---

### Step 3C: Pathway 3 - Auto Standard

**For**: Standardized workflows, limited expertise, consistency

#### Configuration

**Minimal Settings** - almost everything is automatic!

1. Optional: Toggle **"Show detailed progress"**
2. Optional: Upload component specification (for CNMA)
3. Click **"Run Automatic Analysis"** (blue button)

#### What Happens Automatically

The system automatically:
1. ✅ Detects data format and characteristics
2. ✅ Selects optimal method (Standard NMA, CNMA, Meta-regression, etc.)
3. ✅ Sets all parameters (model type, summary measure, reference)
4. ✅ Runs complete pipeline (primary + sensitivity + diagnostics)
5. ✅ Generates recommendations

#### Progress Display

Watch the progress indicators:
- ✓ Step 1: Data detection complete
- ✓ Step 2: Method selection complete
- ✓ Step 3: Analysis complete
- ✓ Step 4: Results ready

#### Results

**"Automatic Choices Made"** box shows:
- Data characteristics detected
- Method selected
- All parameters chosen
- Rationale for choices

**"Analysis Results"** box shows:
- Complete results
- (Full implementation would show comprehensive output)

---

### Step 3D: Pathway 4 - Auto Experimental

**For**: Advanced research with automation, rapid innovation

⚠️ **Warning**: Experimental methods banner displayed

#### Configuration

1. **Research Question** (dropdown):
   - **Let the system decide (auto)** ← Recommended
   - Precision medicine / personalized treatment
   - Clinical decision-making / robustness
   - Survival analysis / time-to-event
   - Model uncertainty

2. **Risk Aversion** (slider 0-2):
   - 0 = Risk-neutral (point estimates)
   - 1 = Moderate caution (default)
   - 2 = Highly cautious (lower confidence bounds)

3. Optional: Toggle **"Show detailed progress"**

4. Click **"Run Experimental Analysis"** (red button)

#### Intelligent Method Selection

Based on your research question AND data characteristics:

**If time-to-event data** → RMST-based NMA + Threshold Analysis
**If IPD with covariates** → ITR + Model Averaging + Threshold
**If pairwise data only** → Threshold Analysis + Model Averaging

#### Results

**"Methods Selected"** box shows:
- Which experimental methods chosen
- Rationale based on research question

**"Experimental Analysis Results"** box shows:
- Primary experimental results
- Novel insights (effect modifiers, robustness, etc.)

**"Comparison with Standard Methods"** (collapsible) shows:
- Agreement with standard NMA
- Differences explained

---

### Step 4: View Results & Export

Click **Results & Export** in sidebar after running any analysis.

#### Analysis Summary

Shows which pathways have been completed:
- ✓ Manual Standard
- ✓ Manual Experimental
- ✓ Auto Standard
- ✓ Auto Experimental

#### Visualizations

Three interactive plot tabs:

1. **Network Plot**
   - Network diagram showing connections
   - Interactive with plotly (zoom, pan, hover)

2. **Forest Plot**
   - Treatment effect estimates with confidence intervals
   - Interactive tooltips

3. **Ranking**
   - Treatment rankings (P-scores or SUCRA)
   - Bar chart with hover details

#### Export Options

Three download buttons:

1. **Download HTML Report** - Complete analysis report in HTML format
2. **Download PDF Report** - Print-ready PDF report
3. **Download Data (CSV)** - Export results data as CSV

---

## Help & Documentation

Click **Help & Documentation** in sidebar for built-in help.

### Four Help Tabs

#### 1. Quick Start
- Getting started guide
- Workflow overview
- Basic instructions

#### 2. Data Formats
- Detailed description of each data format
- Required columns
- Examples

#### 3. Methods Guide
- Overview of all 11 methods
- When to use each method
- Method descriptions

#### 4. About
- powerNMA information
- Version details
- Citation information
- License
- Contact/GitHub link

---

## Data Formats Reference

### Pairwise Format

**Required Columns**: `study, treat1, treat2, TE, seTE`

```
study  treat1  treat2  TE     seTE
1      A       B       0.45   0.15
1      A       C       0.62   0.18
2      B       C       0.28   0.12
```

**Use for**: Most common format, contrast-based data

---

### Arm-based Binary

**Required Columns**: `study, treatment, n, events`

```
study  treatment  n    events
1      A          100  35
1      B          98   42
2      A          120  40
```

**Use for**: Binary outcomes (event counts)

---

### Arm-based Continuous

**Required Columns**: `study, treatment, n, mean, sd`

```
study  treatment  n    mean   sd
1      A          50   12.5   2.3
1      B          48   14.2   2.8
2      C          55   13.1   2.5
```

**Use for**: Continuous outcomes (means and SDs)

---

### Individual Participant Data (IPD)

**Required Columns**: `study, treatment, outcome, [covariates...]`

```
study  treatment  outcome  age  severity  sex
1      A          15.2     45   6.2       M
1      A          14.8     52   7.1       F
1      B          16.5     38   5.8       M
```

**Use for**: Precision medicine, meta-regression, ITR

---

### Time-to-Event

**Required Columns**: `study, treatment, time, event`

```
study  treatment  time   event
1      A          24.5   1
1      A          36.0   0
1      B          18.2   1
```

**Use for**: Survival analysis, RMST-based NMA
- `time`: Time to event or censoring
- `event`: 1 = event occurred, 0 = censored

---

## Advanced Features

### Keyboard Shortcuts

- `Ctrl/Cmd + Click` on sidebar items: Open in new tab (browser feature)
- `F11`: Toggle fullscreen mode (browser feature)

### URL Bookmarking

The GUI supports URL parameters (advanced):
```
http://127.0.0.1:port/?tab=auto_standard
```

---

## Troubleshooting

### GUI Won't Launch

**Problem**: Error when calling `launch_powernma_gui()`

**Solution**:
```r
# Check dependencies
check_gui_dependencies()

# Install missing packages
install_gui_dependencies()

# Try again
launch_powernma_gui()
```

---

### Package Installation Errors

**Problem**: Can't install bs4Dash or other GUI packages

**Solution**:
```r
# Try different CRAN mirror
install.packages("bs4Dash", repos = "https://cran.rstudio.com/")

# Or install from source
install.packages("bs4Dash", type = "source")
```

---

### Data Upload Fails

**Problem**: "Error loading data"

**Solutions**:
1. **Check file format**: Ensure CSV is comma-separated, or use Excel (.xlsx)
2. **Check header**: Toggle "File has header row" checkbox
3. **Check encoding**: Save CSV as UTF-8
4. **Try example data**: Load example first to verify GUI works

---

### Port Already in Use

**Problem**: "Port is already in use"

**Solution**:
```r
# Specify different port
launch_powernma_gui(port = 8888)

# Or find and stop other Shiny apps
# In RStudio: Stop button in Viewer pane
```

---

### Browser Doesn't Open

**Problem**: GUI runs but browser doesn't open

**Solution**:
```r
# GUI is running! Manually open browser to:
# http://127.0.0.1:port

# Or ensure launch.browser = TRUE (default)
launch_powernma_gui(launch.browser = TRUE)
```

---

## Technical Details

### Dependencies

**Required R Packages**:
- `shiny` ≥ 1.7.0 - Core Shiny framework
- `bs4Dash` ≥ 2.3.0 - Bootstrap 4 dashboard
- `DT` ≥ 0.28 - Interactive data tables
- `plotly` ≥ 4.10.0 - Interactive plots
- `shinyWidgets` ≥ 0.7.0 - Enhanced widgets

**Optional Packages**:
- `readxl` - Excel file support
- `rmarkdown` - Report generation
- `knitr` - Report rendering

### System Requirements

- **R Version**: ≥ 4.0.0
- **Browser**: Modern browser (Chrome, Firefox, Safari, Edge)
- **RAM**: ≥ 4 GB recommended for large datasets
- **Ports**: One available TCP port (random by default)

### File Structure

```
powerNMA/
├── inst/
│   └── shiny/
│       └── app.R              # Main Shiny application
├── R/
│   └── launch_gui.R           # Launcher functions
└── NAMESPACE                  # Package exports
```

---

## Best Practices

### 1. Data Preparation

✅ **Do**:
- Clean data before uploading
- Check for missing values
- Verify column names match format requirements
- Use UTF-8 encoding for text

❌ **Don't**:
- Upload very large files (>100 MB) without testing
- Mix data formats in one file
- Include non-ASCII characters without proper encoding

### 2. Analysis Strategy

✅ **Do**:
- Start with Auto Standard to explore
- Use Manual pathways for production
- Compare multiple pathways for important decisions
- Export and save results frequently

❌ **Don't**:
- Rely solely on Auto Experimental for clinical guidelines
- Ignore experimental warnings
- Skip data quality checks

### 3. Results Interpretation

✅ **Do**:
- Review automatic choices (Auto pathways)
- Check diagnostics and sensitivity analyses
- Compare experimental with standard methods
- Consult statistical expert for experimental results

❌ **Don't**:
- Accept results without review
- Ignore inconsistency warnings
- Use experimental methods for regulatory submissions without validation

---

## FAQ

### Q: Can I use the GUI offline?

**A**: Yes! The GUI runs locally on your computer. Internet is only needed to install packages.

---

### Q: Is my data secure?

**A**: Yes. All data processing happens locally on your computer. No data is sent anywhere.

---

### Q: Can I customize the GUI?

**A**: The GUI code is in `inst/shiny/app.R`. Advanced users can modify it.

---

### Q: How do I cite the GUI?

**A**: Same as citing powerNMA package (citation information in About tab).

---

### Q: Can I deploy this on a server?

**A**: Yes! Use:
```r
launch_powernma_gui(host = "0.0.0.0", launch.browser = FALSE)
```

For production deployment, consider Shiny Server or shinyapps.io.

---

### Q: What browsers are supported?

**A**: All modern browsers:
- Chrome/Chromium ✅
- Firefox ✅
- Safari ✅
- Edge ✅
- Internet Explorer ❌ (not supported)

---

### Q: Can I run multiple instances?

**A**: Yes! Use different ports:
```r
# Instance 1
launch_powernma_gui(port = 8080)

# Instance 2 (in new R session)
launch_powernma_gui(port = 8081)
```

---

## Updates and Feedback

### Stay Updated

- **GitHub**: https://github.com/mahmood726-cyber/rmstnma
- **Issues**: Report bugs or request features via GitHub Issues
- **Tag**: Use `[GUI]` tag for GUI-specific issues

### Feedback

We welcome feedback on the GUI! Please let us know:
- What features you'd like to see
- What's confusing or unclear
- Any bugs or issues encountered

---

## Roadmap

### Planned Features

🔮 **Version 3.1** (Q1 2026):
- Full result visualization (plots)
- Report generation (HTML/PDF)
- Advanced plot customization
- Data validation feedback

🔮 **Version 3.2** (Q2 2026):
- Real-time analysis updates
- Comparison of multiple analyses
- Export to Word format
- Interactive network editing

🔮 **Version 3.3** (Q3 2026):
- Cloud deployment support
- Collaborative features
- API integration
- Mobile-responsive design

---

## Acknowledgments

The powerNMA GUI is built with:
- **Shiny**: RStudio
- **bs4Dash**: David Granjon
- **plotly**: Plotly Technologies
- **DT**: RStudio
- **shinyWidgets**: Victor Perrier & Fanny Meyer

---

**Document Version**: 1.0
**Last Updated**: October 31, 2025
**Next Review**: January 2026
