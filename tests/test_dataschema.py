import pytest
from strictyaml.yamllocation import YAMLChunk as YC
from wingstructure.data import data_schema as ds


def tostr(adict):
    resdict = {}
    for key, val in adict.items():
        resdict[key] = str(val)
    return resdict


class TestGeometry():
    """
    Tests geometry schemas
    """

    def test_vector(self):
        oneval = {'x': 10.0}
        filled = {'x': 10.0, 'y': 0.0, 'z': 0.0}
        full = {'x': 0.5, 'y': 10.1, 'z': -0.4}

        def vector_validate(val):
            return dict(ds.vector.validate(YC(tostr(val))))

        assert vector_validate(oneval) == filled
        assert vector_validate(full) == full

    def test_control_surface(self):

        def cs_validate(val): dict(ds.control_surface.validate(YC(tostr(val))))

        valid_cs = {'name': 'aileron1',
                    'type': 'aileron',
                    'span-start': 0.5,
                    'span-end': 1.0}
        assert cs_validate(valid_cs) == valid_cs
