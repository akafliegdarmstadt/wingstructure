wingdict = {
    'pos': {'x': 0.0, 'y':0.0, 'z':0.0},
    'rot': {'x': 0.0, 'y':0.0, 'z':0.0},
    'sections':[
        {
            'pos': {'x': 0.0, 'y':0.0, 'z':0.0},
            'twist': 0.0,
            'chord': 0.9,
            'airfoil': 'Clark-y'
        },
        {
            'pos': {'x': 0.0, 'y':1.0, 'z':0.0},
            'twist': 0.0,
            'chord': 0.7,
            'airfoil': 'Clark-y'
        },
        {
            'pos': {'x': 0.0, 'y':5.0, 'z':0.0},
            'twist': -0.1,
            'chord': 0.3,
            'airfoil': 'Clark-y'
        }
    ],
    'control-surfaces':[
        {
            'name':'spoiler1',
            'type':'spoiler',
            'span-start': 0.7,
            'span-end': 1.4,
        },
        {
            'name':'aileron1',
            'type':'aileron',
            'span-start': 0.7,
            'span-end': 1.4,
            'chord-pos': 0.7
        }
    ]
}