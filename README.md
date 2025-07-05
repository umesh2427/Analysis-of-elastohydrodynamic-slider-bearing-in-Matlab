# Analysis-of-elastohydrodynamic-slider-bearing-in-Matlab
Using Reynolds and energy equations, slider bearing is analyzed with the help of the finite difference method in Matlab to obtain performance parameters such as pressure and temperature.
# Elastohydrodynamic Analysis of Slider Bearing in MATLAB 🛠️

A comprehensive MATLAB simulation of an **elastohydrodynamic slider bearing**, showcasing fluid–structure interactions under varying load conditions using the Reynolds equation and structural elasticity coupling.

## 📘 Overview

This project models a slider bearing where the lubricant film behaves elastohydrodynamically—taking into account both **hydrodynamic pressure** and **deformation** of the bearing surfaces. The MATLAB code captures how varying loads, bearing stiffness, and temperature affect film thickness, pressure distribution, and load-carrying capacity, based on established academic formulations :contentReference[oaicite:1]{index=1}.

## 🚀 Features

- **Reynolds equation** solver for isothermal hydrodynamic film pressure
- Elastic deformation analysis of bearing surfaces via coupled elasticity equations
- Parametric study across:
  - Applied load / bearing width
  - Elastic moduli (rigid vs. flexible housing)
  - Oil film thickness and pressure distribution
- Visualization tools:
  - Oil film thickness maps
  - Pressure contour plots
  - Film thickness and bearing pressure vs. load graphs

## ⚙️ Inputs & Configuration

- **Geometric parameters**: bearing length, width, clearance  
- **Material properties**: Young’s modulus, Poisson’s ratio  
- **Lubricant properties**: viscosity, temperature  
- **Load conditions**: range of applied loads, surface stiffness variations



## 🧠 Key Results

- Deformation leads to **non-uniform pressure distribution**, often skewed toward one edge.  
- The **minimum film thickness decreases** as load increases — risk of mixed or boundary lubrication regimes.  
- Flexible bearing housings show **edge-localized thin films**, while rigid ones distribute loads more evenly, echoing literature findings :contentReference[oaicite:2]{index=2}.

## 🎯 Applications & Significance

- Ideal for **tribology studies** and understanding film lubrication.  
- Useful for **designers** of bearings in engines, gearboxes, and heavy machinery.  
- Educational tool for illustrating **fluid–structure coupling** in mechanical systems.

## 🔮 Future Enhancements

- Model **thermal effects** (temperature-pressure dependent viscosity).  
- Introduce **transient or dynamic loads** to analyze mixed lubrication.  
- Integrate **surface deformation via FEM** or advanced elasticity models.  
- Expand study to other bearing types (journal, tilting-pad).



---


