import collections
import strictyaml as sy

# helper definitions
vector = sy.Map({
    sy.Optional('x', default=0.0): sy.Float(),
    sy.Optional('y', default=0.0): sy.Float(),
    sy.Optional('z', default=0.0): sy.Float()
})

zeros = collections.OrderedDict({'x': 0.0, 'y': 0.0, 'z': 0.0})

scaling = sy.OrValidator(vector, sy.Float())

# geometry definitions

wingsection = sy.Map({
    'pos': vector,
    sy.Optional('twist', default=0.0): sy.Float(),
    'chord': sy.Float(),
    'airfoil': sy.Str()
})

control_surface = sy.Map({
    'type': sy.Enum(['aileron',
                  'flap',
                  'flaperon',
                  'spoiler',
                  'elevator',
                  'rudder']),
    'span-start': sy.Float(),
    'span-end': sy.Float(),
    sy.Optional('chord-pos'): sy.OrValidator(sy.Float(),
                 sy.FixedSeq([sy.Float(), sy.Float()]))
})

wing = sy.Map({
    sy.Optional('pos', default=zeros): vector,
    sy.Optional('rot', default=zeros): vector,
    sy.Optional('scale', default=1.0): scaling,
    'sections': sy.Seq(wingsection),
    sy.Optional('control-surfaces'): sy.MapPattern(sy.Str(), control_surface)
})

geometry = sy.Map({'wing': wing,
                   sy.Optional('elevator'): wing,
                   sy.Optional('fuselage'): sy.Any()})

# aerodynamic definition



# overall definition

data = sy.Map({
    'geometry': geometry,
    sy.Optional('aerodynamic'): sy.Any()
    })


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
    
    loaded = sy.dirty_load(yamlstr, data, allow_flow_style=True)
    
    return clear(loaded)


def clear(data):
    '''helper function removing metadata'''

    if type(data) == sy.representation.YAML:
        return clear(data.data)
    elif type(data) == collections.OrderedDict:
        return {key: clear(value) for key, value in data.items()}
    elif type(data) == list:
        return [clear(item) for item in data]
    else:
        return data


if __name__ == '__main__':
    loaded = loaddata('../../tests/test.yaml')
    print(loaded.__repr__())
