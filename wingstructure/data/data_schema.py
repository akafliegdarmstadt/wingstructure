from collections import OrderedDict
from strictyaml import load, Map, Str, Int, Float, Seq, FixedSeq,\
                        Any, Enum, Optional, OrValidator, YAMLError,\
                        representation

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
                Optional('elevator'): wing})

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
    
    return clear(loaded)


def clear(data):
    '''helper function removing metadata'''

    if type(data) == representation.YAML:
        return clear(data.data)
    elif type(data) == OrderedDict:
        return {key: clear(value) for key, value in data.items()}
    elif type(data) == list:
        return [clear(item) for item in data]
    else:
        return data


if __name__ == '__main__':
    loaded = loaddata('../../tests/test.yaml')
    print(loaded.__repr__())
