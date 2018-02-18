# wingstructure
a python tool for wing structure calculations

![lift distribution](examples/Liftdistribution.png)

## Usage
Use pip to install the package:
```sh
pip install https://github.com/helo9/wingstructure/archive/master.zip
```

### Python
Import relevant classes:
```python
from wingstructure import Wing, LiftAnalysis, LiftAndMomentAnalysis
```

Create simple geometry:
```python
span_positions = [0, 2, 5, 7]
chord_lengths = [1, 0.9, 0.6, 0.3]
offsets = [0, 0.1, 0.4, 0.7]
twists = [0]*4
airfoils = [None]*4
wing = Wing.create_from_planform(span_positions, chord_lengths, offsets, twists, airfoils)
wing.set_root_pos(0.0)
wing.set_airbrake(1.5,2.9)
```

![geometry](examples/wing.png)

Create Analysis
```python
liftana = LiftAnalysis(wing)
```

Calculate Distribution of lift Coefficient
```python
α, distribution = liftana.calculate(lift=0.8)
```
