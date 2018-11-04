from collections import OrderedDict
from strictyaml import load, Map, Str, Int, Float, Seq, FixedSeq,\
                        Any, Enum, Optional, OrValidator, YAMLError

# helper definitions
vector = Map({
    Optional('x', default=0.0): Float(),
    Optional('y', default=0.0): Float(),
    Optional('z', default=0.0): Float()
})

zeros = OrderedDict({'x': 0.0, 'y': 0.0, 'z': 0.0})

scaling = OrValidator(vector, Float())


# geometry definitions

wingsection = Map({
    'pos': vector,
    Optional('twist', default=0.0): Float(),
    'chord': Float(),
    'airfoil': Str()
})

control_surface = Map({
    'name': Str(),
    'type': Enum(['aileron',
                  'flap',
                  'flaperon',
                  'spoiler',
                  'elevator',
                  'rudder']),
    'span-start': Float(),
    'span-end': Float(),
    'chord-pos': OrValidator(Float(),
                             FixedSeq([Float(), Float()]))
})

wing = Map({
    Optional('pos', default=zeros): vector,
    Optional('rot', default=zeros): vector,
    Optional('scale', default=1.0): scaling,
    'sections': Seq(wingsection),
    Optional('control-surfaces'): Seq(control_surface)
})

geometry = Map({'wing': wing,
                'elevator': wing})

# overall definition

data = Map({'geometry': geometry})


def loaddata(filename: str):
    """opens file and trys reading construction data
    
    Parameters
    ----------
    filename : str
        file containing construction data
    
    Returns
    -------
    OrderedDict
        read data
    """

    with open(filename, 'r') as yamlfile:
        yamlstr = yamlfile.read()
    
    loaded = load(yamlstr, data)
    
    return loaded


if __name__ == '__main__':
    loaded = loaddata('../../tests/test.yaml')
    print(loaded.__repr__())
