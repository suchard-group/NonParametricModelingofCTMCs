# Non-Parametric Modeling of Continuous-Time Markov Chains

[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL_v2.1-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
[![BEAST](https://img.shields.io/badge/BEAST-hmc--clock-green)](https://github.com/beast-dev/beast-mcmc)
[![BEAGLE](https://img.shields.io/badge/BEAGLE-hmc--clock-orange)](https://github.com/beagle-dev/beagle-lib)

This repository contains the scripts and XML files required to reproduce the analyses presented in the paper  
**“NonParametric Modeling of Continuous-Time Markov Chains”**  
by *Filippo Monti*, *Xiang Li*, *Marc A. Suchard*.

---

## 📂 Folder Descriptions
### `code/xml/`
This folder contains the XML files used to perform the analyses described in the manuscript.  
Each XML file produces a log file that stores the MCMC samples for the parameters of interest.  
The XML files are executed using **BEAST** and **BEAGLE**, whose installation is described below.

There are **two XML files per dataset**, one using the log-linear model and one using the Gaussian process model.  
The three datasets analyzed are:
- A dataset **simulated** using [PiBuss](https://beast.community/pibuss) (available directly within BEAST)
- A dataset of **rabies virus sequences** sampled from North American bats
- A dataset of **influenza H3N2 sequences** sampled from the United States  

For the influenza dataset, an additional file contains a set of **empirical trees** sampled from the posterior distribution inferred in a previous analysis as well as an XML to generate a **world map**.

### `code/R_code/`
This folder contains the R scripts used to analyze the log files and generate the figures presented in the manuscript.

#### Contents
- **`NonParametricModelingofCTMCs.Rproj`**  — R project file for easy loading of the repository in RStudio. If the files are opened and run within this project, all paths will be adapted automatically.

- **`R_figures/`** — Four Quarto files (`.qmd`) that generate all manuscript figures.  
  There is one Quarto file per dataset and one that reproduces the figures on gradient evaluation time comparison. 
  In the first code chunk the Quarto files *automatically* load all other scripts in the `R_code` folder and its subfolders, provided the folder hierarchy is maintained.
  The first section named **Preliminaries** contains two code chunks: one to set the file names and the dummy variable **save_plots** (if TRUE, the figures will be saved to the output folder), and one that loads the required fonts (if getting an error here, check the file `/code/R_code/R_functions/loadFonts.R`).

- **`path.R`** — Defines the working directories. The paths are not computer-specific, so they should not be changed unless a different (within project) folder hierarchy is desired.
- **`libraries.R`** — Loads all required R packages.

- **`R_classes/`** and **`R_functions/`** - Contain the R classes and functions used by the Quarto documents. These are sourced automatically if the directory structure is preserved.

---

## ⚙️ Requirements

To reproduce the analyses, you will need:
- **BEAST (hmc-clock branch)**
- **BEAGLE (v4.0.0 or hmc-clock branch)**
- **R** (≥4.0.0)
- All R packages listed in `R_code/libraries.R`

The R packages are automatically loaded by the Quarto files in `R_code/R_figures`.

---

## 🧩 Installing BEAST and BEAGLE

### BEAGLE

You may use either:
- The v4.0.0 release: https://github.com/beagle-dev/beagle-lib/releases/tag/v4.0.0  
- Or the `hmc-clock` branch: https://github.com/beagle-dev/beagle-lib

#### macOS

    xcode-select --install
    brew install libtool autoconf automake
    git clone https://github.com/beagle-dev/beagle-lib.git
    cd beagle-lib
    git checkout hmc-clock
    mkdir build && cd build
    cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF ..
    sudo make install

#### Linux

    sudo apt-get install build-essential autoconf automake libtool git pkg-config openjdk-9-jdk
    git clone https://github.com/beagle-dev/beagle-lib.git
    cd beagle-lib
    git checkout hmc-clock
    mkdir build && cd build
    cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF ..
    sudo make install

> 🗂️ Libraries are typically installed into `/usr/local/lib`.

#### Setting up BEAGLE
If opting not to use the [v4.0.0 release of BEAGLE](https://github.com/beagle-dev/beagle-lib/releases/tag/v4.0.0), please follow the [BEAGLE installation instructions](https://github.com/beagle-dev/beagle-lib), but be sure to get the `hmc-clock` branch.

For Mac users, the following commands will compile the CPU version of BEAGLE.
Follow the [instructions](https://github.com/beagle-dev/beagle-lib) if you need to install any other dependent software; ignore the first 2 lines if you already have all requisite dependencies installed.

```
xcode-select --install
brew install libtool autoconf automake
git clone https://github.com/beagle-dev/beagle-lib.git
cd beagle-lib
git checkout hmc-clock
mkdir build
cd build
cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF ..
sudo make install
```

For Linux users, the commands are similar.

```
sudo apt-get install build-essential autoconf automake libtool git pkg-config openjdk-9-jdk
git clone https://github.com/beagle-dev/beagle-lib.git
cd beagle-lib
git checkout hmc-clock
mkdir build
cd build
cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF ..
sudo make install
```

The libraries are installed into `/usr/local/lib`.

### BEAST

#### Compilation

    git clone https://github.com/beast-dev/beast-mcmc.git
    cd beast-mcmc
    git checkout hmc-clock
    ant

#### Installing `ant`

- macOS: `brew install ant` using [Homebrew](https://brew.sh/)
- Linux: `sudo apt-get install ant`

After compilation, the file `beast.jar` will be located at:

    beast-mcmc/build/dist/beast.jar

---

## ▶️ Running BEAST with BEAGLE

To verify BEAGLE installation, run:

    java -jar /path/to/beast-mcmc/build/dist/beast.jar -beagle_info

If BEAST cannot locate BEAGLE, try:

    java -Djava.library.path=/path/to/beagle \
      -jar /path/to/beast-mcmc/build/dist/beast.jar -beagle_info

Or add BEAGLE to your library path:

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

More information is available at: https://github.com/beagle-dev/beagle-lib

---
## 🧪 Reproducing the Analyses

### 🧬 Simulation

Generate log files for the **log-linear** and **Gaussian process** models.

**Log-Linear Model**

    java -Djava.library.path=/usr/local/lib \
      -jar /path/to/beast-mcmc/build/dist/beast.jar \
      -seed 666 -overwrite \
      ../xmls/simulation_LL.xml

**Gaussian Process Model**

    java -Djava.library.path=/usr/local/lib \
      -jar /path/to/beast-mcmc/build/dist/beast.jar \
      -seed 666 -overwrite \
      ../xmls/simulation_GP.xml

---

### 🦇 Rabies Virus in North American Bats

**Log-Linear Model**

    java -Djava.library.path=/usr/local/lib \
      -jar /path/to/beast-mcmc/build/dist/beast.jar \
      -seed 666 -overwrite \
      ../xmls/rabies_LL.xml

**Gaussian Process Model**

    java -Djava.library.path=/usr/local/lib \
      -jar /path/to/beast-mcmc/build/dist/beast.jar \
      -seed 666 -overwrite \
      ../xmls/rabies_GP.xml

---

### 🧫 Influenza H3N2

This analysis requires a set of **empirical trees** `code/xml/influenza_empirical_trees`.

[//]: # (Extract them into the `xml` folder using:)

[//]: # ()
[//]: # (    tar -xzf ../xmls/airCommunitiesMM_500.trees.tar.gz -C ../xmls/)

[//]: # (Then run:)

**Log-Linear Model**

    java -Djava.library.path=/usr/local/lib \
      -jar /path/to/beast-mcmc/build/dist/beast.jar \
      -seed 666 -overwrite \
      ../xmls/influenza_LL.xml

**Gaussian Process Model**

    java -Djava.library.path=/usr/local/lib \
      -jar /path/to/beast-mcmc/build/dist/beast.jar \
      -seed 666 -overwrite \
      ../xmls/influenza_GP.xml

---

### 🌍 Generating MCC trees 
To generate MCC trees with **TreeAnnotator** (included with BEAST) use:

    treeannotator -heights mean \
      -burnin 10 \
      input.trees \
      output.tree

where `input.trees` is the file containing the posterior trees generated by BEAST and `output.tree` is the name of the output MCC tree file. The R files assumes the MCC trees are named `rabies_GP_MCC.tree` and `influenza_GP_MCC.tree` and are stored in the `./output` directory.

## 📖 Citation

If you use this repository, please cite:

> Monti, F., Li, X., and Suchard, M.A.  
> *NonParametric Modeling of Continuous-Time Markov Chains.*  
> (Manuscript submitted).

---

## 🧠 Authors

- **Filippo Monti** — UCLA  
- **Xiang Li** — Tulane University
- **Marc A. Suchard** — UCLA

---

## 🪪 License

This repository is released under the LGPL v2.1 License:  
https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html

