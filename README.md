# Theo Jansen Simulation

This note explains the Python file `theo_jansen.py`.
It follows the mechanism topology documented in:

- Shin, Deshpande, and Sulzer, "Design of a Single Degree-of-Freedom, Adaptable Electromechanical Gait Trainer for People With Neurological Injury", *Journal of Mechanisms and Robotics*, 2018.
- Jadav et al., "Kinematic Performance of a Customizable Single Degree-of-Freedom Gait Trainer", supplied as `sn-article.pdf`.

The experimental link-length table from the later paper is:

| link       |  L1 |  L2 |  L3 |  L4 |   L5 |   L6 |   L7 |   L8 |  L9 | L10 |  L11 |  L12 |
| ---------- | --: | --: | --: | --: | ---: | ---: | ---: | ---: | --: | --: | ---: | ---: |
| nominal cm |  11 |  45 |  36 |  34 | 48.5 | 41.5 | 60.5 | 41.5 |  42 |  43 | 26.5 | 54.5 |

The script uses the paper's example adjustable values:

`L1 = 11.29 cm`, `L4 = 32.93 cm`, `L8 = 41.78 cm`.

With the gait/world frame rotated by `-11.4 deg` relative to the schematic ground-link frame, the script produces:

`x-span = 50.08 cm`, `y-span = 12.74 cm`.

This matches the paper-reported simulated output for the chosen human gait envelope, approximately `[50.10, 12.75] cm`.

---

## 1. Coordinate Frames

The mechanism is solved first in a mechanism frame `M`.

- `P0 = [0, 0]^T` is the crank ground pivot.
- `P3 = [L4, 0]^T` is the second ground pivot.
- Link `L4` is the fixed ground link.
- Link `L1` is the input crank.
- The input angle `theta1` is swept from `0` to `2*pi`.

The gait/world frame `G` is obtained by rotating all solved points by:

`alpha = -11.4 deg`.

That rotation does not change the mechanism closure. It only expresses the ankle path in the same horizontal/vertical gait axes used for stride length and step height.

---

## 2. Link Direction Convention

Define:

`e(theta) = [cos(theta), sin(theta)]^T`.

The script uses link angles `theta1 ... theta12` with these vector directions:

- `L1`: `P0 -> P1`
- `L2`: `P1 -> P2`
- `L3`: `P3 -> P2`
- `L4`: `P0 -> P3`, fixed, so `theta4 = 0`
- `L5`: `P2 -> P4`
- `L6`: `P3 -> P4`
- `L7`: `P1 -> P5`
- `L8`: `P3 -> P5`
- `L9`: `P5 -> P6`
- `L10`: `P4 -> P6`
- `L11`: `P5 -> PE`
- `L12`: `PE -> P6`

The end-effector / ankle point is `PE`, the joint between `L11` and `L12`.

The unknown angle vector at each crank angle is:

`q = [theta2, theta3, theta5, theta6, theta7, theta8, theta9, theta10, theta11, theta12]^T`.

`theta1` is prescribed. `theta4` is fixed.

---

## 3. Closed Vector Loops

The mechanism has five closure equations, each with x and y components, giving 10 scalar equations for 10 unknown angles.

### Upper four-bar: L1-L2-L3-L4

Path: `P0 -> P1 -> P2 -> P3 -> P0`

`L1*e1 + L2*e2 - L3*e3 - L4*e4 = 0`

### Lower four-bar: L1-L7-L8-L4

Path: `P0 -> P1 -> P5 -> P3 -> P0`

`L1*e1 + L7*e7 - L8*e8 - L4*e4 = 0`

### Coupler triangle: L3-L5-L6

Path: `P3 -> P2 -> P4 -> P3`

`L3*e3 + L5*e5 - L6*e6 = 0`

### Parallelogram-like loop: L6-L8-L9-L10

Path: `P3 -> P5 -> P6 -> P4 -> P3`

`L8*e8 + L9*e9 - L6*e6 - L10*e10 = 0`

The paper calls this a parallelogram mechanism. The optimized lengths are close but not exactly equal in opposite pairs, so the code treats it as a general four-bar closure.

### Foot triangle: L9-L11-L12

Path: `P5 -> PE -> P6 -> P5`

`L11*e11 + L12*e12 - L9*e9 = 0`

Together:

`F(q, theta1, L) = 0`,

implemented in the function `loop_equations`.

---

## 4. Forward Kinematics

Once the angles are solved, joint positions are computed as:

```python
P0 = np.array([0.0, 0.0])
P3 = np.array([L4,  0.0])
P1 = P0 + L1 * e(theta1)
P2 = P1 + L2 * e(theta2)
P5 = P1 + L7 * e(theta7)
P4 = P2 + L5 * e(theta5)
P6 = P5 + L9 * e(theta9)
PE = P5 + L11 * e(theta11)   # ankle / end-effector
```

The equivalent alternative definitions,

```python
P2 = P3 + L3 * e(theta3)
P5 = P3 + L8 * e(theta8)
P4 = P3 + L6 * e(theta6)
P6 = P4 + L10 * e(theta10)
PE = P6 - L12 * e(theta12)
```

are enforced by the vector loops and serve as closure residual checks.

---

## 5. Numerical Solving

For each `theta1`:

1. A geometric circle-intersection assembly is used only as a branch seed for the upper four-bar.
2. The full nonlinear loop system is solved using `scipy.optimize.fsolve`.
3. The previous solution is used as the initial guess for the next step, which preserves the physical assembly mode through a full revolution.

A damped Newton solver with finite-difference Jacobians is included as a fallback. It solves the same equations and achieves residuals near `1e-14 cm` during verification.

---

## 6. Plots Produced

Running the script generates:

- A composite figure similar to Fig. 4 of the paper: mechanism sketch, endpoint markers, and x/y gait-cycle curves.
- A 3x3 validation grid with reference-vs-simulation endpoint paths and RMSE labels.
- End-effector trajectory `x` versus `y`.
- `x` and `y` versus gait cycle percent.
- End-effector velocity curves for constant crank speed.
- Closed-loop residual validation across the full crank revolution.
- Mechanism snapshots at representative crank angles.
- Sensitivity of stride length and step height to adjustable links `L1`, `L4`, and `L8`.
- Least-squares span-to-link mapping demonstration.
- A live mechanism animation when run in an interactive environment.

---

## 7. Effect of Adjustable Links

**Link `L1` (crank radius):** Increasing it generally increases excitation of both four-bars, changing both stride length and vertical lift. It is the primary amplitude control for the gait envelope.

**Link `L4` (ground pivot separation):** Changing it shifts the base geometry of both four-bars simultaneously, strongly affecting horizontal span and the timing and shape of the swing arc. It is the dominant parameter for stride length tuning.

**Link `L8` (lower four-bar follower):** Belongs to both the lower four-bar and the parallelogram-like loop. It strongly relocates joint `P5`, which then moves the entire foot triangle. It is the principal control for step height and end-effector loop shape.

This is why the papers focus on `L4` and `L8`, and the later paper additionally makes `L1` adjustable. Together these three links provide a compact, clinically meaningful parameterization of the gait envelope.

---

## 8. Least-Squares Span Mapping

The paper describes an empirical map from desired gait spans to adjustable link lengths.

Let:

- `Lambda in R^(3 x n)` store sampled `[L1; L4; L8]` vectors.
- `Sigma in R^(2 x n)` store the corresponding simulated `[xspan; yspan]` vectors.
- `Psi in R^(3 x 2)` be the linear map from spans to link lengths.

The least-squares solution is:

`Psi = Lambda * pinv(Sigma)`.

For a new desired span:

`[L1; L4; L8] ~= Psi * [xspan; yspan]`.

This is implemented in the function `run_least_squares_span_map`. The result gives a direct starting estimate that can be refined by a nonlinear polishing step if sub-centimeter accuracy is required.

---

## 9. Important Accuracy Note

The source PDFs do not include the raw human ankle marker data or the authors' original code.
The blue reference curve included in validation plots is therefore a smooth gait-like envelope reconstructed for visual comparison only. The reproduced research contribution is the closed-chain mechanism kinematics and the output gait envelope. With the supplied example adjustable lengths, the simulated stride and step height spans match the paper-reported values to within `0.02 cm`.
