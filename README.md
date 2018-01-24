# wingstructure
python tool for wing structure calculations

![lift distribution](examples/Liftdistribution.png)

## Usage
Use pip to install Package
```sh
pip install https://github.com/helo9/wingstructure/archive/master.zip
```

### Python
Import relevant classes
```python
from wingstructure import WingExt, LiftAnalysis, LiftAndMomentAnalysis
```

Create simple geometry
```python
span_positions = [0, 2, 5, 7]
chord_lengths = [1, 0.9, 0.6, 0.3]
offsets = [0, 0.1, 0.4, 0.7]
twists = [0]*4
airfoils = [None]*4
wing = WingExt.create_from_planform(span_positions, chord_lengths, offsets, twists, airfoils)
wing.set_root_pos(0.0)
wing.set_airbrake(1.5,2.9)
wing.set_flap('flap', 2, 5,[0.3,0.3])
wing.set_flap('flap2', 5, 7, [0.3,0.2])
```

![geometry](examples/wing.png)

Create Analysis
```python
liftana = LiftAnalysis(wing)
```

Calculate Distribution for Lift Coefficient
```python
α, distribution = liftana.calculate(lift=0.8)
```
