# 🦿 Simulation of a 12-Link Jansen Walking Mechanism

## 📌 Overview

This project focuses on modeling and simulating a **12-bar Theo Jansen linkage** that produces a walking-like motion. The mechanism operates with a **single rotational input**, making it efficient for applications such as assistive gait systems and robotic walkers.

---

## 🎯 Purpose

* Construct a kinematic representation of a multi-link system
* Solve constraint equations for closed-loop mechanisms
* Generate and analyze a walking trajectory
* Evaluate stride length and lift height
* Cross-check results with published studies

---

## 📚 Source Material

The formulation is based on:

* Shin et al. (2018): Adaptive electromechanical gait system
* Jadav et al.: Performance study of customizable walking linkages

---

## ⚙️ Mechanism Details

The structure includes:

* 12 rigid links connected through joints
* 5 independent closed kinematic loops
* One driving crank (input)

### Important Points:

* Fixed joints: P0 and P3
* Output node: PE (foot position)
* Input variable: θ₁ (crank rotation)

---

## 📐 Governing Equation

The motion of the system is defined by:

```id="eqA"
F(q, θ₁, L) = 0
```

Where:

* `q` represents unknown angular positions
* `θ₁` is the known input angle
* `L` denotes link lengths

---

## 🔁 Constraint Equations

### First Loop (Upper Chain)

```id="eqB"
L1e1 + L2e2 − L3e3 − L4e4 = 0
```

### Second Loop (Lower Chain)

```id="eqC"
L1e1 + L7e7 − L8e8 − L4e4 = 0
```

### Third Loop (Triangular Section)

```id="eqD"
L3e3 + L5e5 − L6e6 = 0
```

### Fourth Loop (Coupled Linkage)

```id="eqE"
L8e8 + L9e9 − L6e6 − L10e10 = 0
```

### Fifth Loop (Foot Assembly)

```id="eqF"
L11e11 + L12e12 − L9e9 = 0
```

These equations ensure that all links remain connected and consistent during motion.

---

## 🧮 Position Computation

Joint coordinates are obtained step-by-step:

```id="eqG"
P1 = P0 + L1e(θ₁)
P2 = P1 + L2e(θ₂)
P5 = P1 + L7e(θ₇)
PE = P5 + L11e(θ₁₁)
```

Other equivalent formulations are satisfied through loop constraints.

---

## 🛠️ Solution Strategy

To compute unknown angles:

* Use MATLAB’s `fsolve` where available
* Otherwise apply an iterative Newton-based method
* Use the previous step’s solution as the initial guess
* Maintain consistent motion branch across full rotation

---

## 📊 Outputs Generated

The simulation provides:

* Foot path in Cartesian space
* Motion variation over gait cycle
* Velocity profiles
* Animated linkage movement
* Error verification of constraints
* Effect of parameter variation

---

## 📏 Observed Results

Typical simulation output:

* **Step length ≈ 50 cm**
* **Step height ≈ 12.7 cm**

These values closely match those reported in literature.

---

## 🔧 Role of Key Links

| Link | Contribution                              |
| ---- | ----------------------------------------- |
| L1   | Controls overall motion input             |
| L4   | Defines base geometry and stride          |
| L8   | Influences vertical motion and trajectory |

---

## 📈 Mapping Desired Motion to Link Lengths

To estimate link dimensions from required motion:

```id="eqH"
Ψ = Λ · pinv(Σ)
```

This provides an approximate relationship between gait parameters and link configuration.


---
## 👥 Team Members
- Rutva Mehta - 202401116
- Yash Ravat - 202401182 
- Om Pandya - 202401136  
- Harsh Asnani - 202401062  
- Het Lathiya - 202401104

---
## 📁 Drive Folder
This Drive folder includes all the source code used in the demonstration of our project.
Link - https://drive.google.com/drive/folders/1SgetUaWjeaJx1NPeDUs9qbj-JOlweufU


---

