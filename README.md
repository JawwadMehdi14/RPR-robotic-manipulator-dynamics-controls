# 🤖 3-DoF Fixed-Base Robot Manipulator — Kinematics, Dynamics & Control

## 📌 Overview
This project presents the complete modeling, analysis, and control of a **fixed-base 3 Degrees of Freedom (DoF) robot manipulator**, featuring:
- **Two prismatic joints** (linear motion)
- **One revolute joint** (rotational motion)

Using **MATLAB**, **Simulink**, and the **Robotics Toolbox**, the robot was analyzed for:
- **Kinematics** (forward & inverse)
- **Differential kinematics**
- **Dynamics** (via Lagrange and Newton–Euler methods)
- **Control strategies** (joint space & operational space)

The work integrates theoretical modeling with practical simulation to evaluate and compare various controllers under trajectory tracking and force control tasks.

---

## 🎯 Objectives
1. **Model the Robot** using the Denavit–Hartenberg (DH) convention.
2. **Analyze Kinematics**:
   - Forward & inverse kinematics using transformation matrices and symbolic equations.
   - Verification using the MATLAB Robotics Toolbox.
3. **Perform Differential Kinematics**:
   - End-effector velocity analysis via geometric and analytical Jacobians.
4. **Derive Dynamics**:
   - Using both Lagrange and Newton–Euler formulations for accuracy comparison.
5. **Implement Controllers**:
   - PD Control
   - Inverse Dynamics Control
   - Adaptive Control
   - Force Control
   - Compliance Control
   - Impedance Control
   - Admittance Control
   - Parallel Force/Position Control
6. **Evaluate Performance** in both **joint space** and **operational space** through MATLAB Simulink simulations.

---

## 🛠️ Features
- **Comprehensive Robot Modeling** using the DH convention.
- **Symbolic and Numerical Analysis** of kinematics and dynamics.
- **Multiple Control Strategies**:
  - Joint space controllers.
  - Operational space controllers.
  - Direct and indirect force control.
- **Disturbance Testing** for robustness analysis.
- **Simulation-based Validation** using MATLAB Simulink.
- **Controller Performance Comparison** under various scenarios.

---
