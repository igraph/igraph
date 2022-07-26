# This file is used to test the handling of decimal points
# under different locales.
graph [
    node [
        id 1
        value 1.23
    ]
        node [
        id 2
        value 4.87e-5
    ]
    edge [
        source 1
        target 2
        weight -3.779e+2
    ]
]
