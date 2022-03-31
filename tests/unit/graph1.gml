# comment lines are ignored
Version 1
graph [
    directed 1
    node [
        id 1
        a [ b 1 ]
        b 2
        c 5.6
        graphics [ x 0 y 0 ]
    ]
    node [
        a 1
        id 2
        a 2
        b "asd"
        d -Inf
    ]
    edge [ 
        source 5 target 1 
        weight -0.45
        label -0.45
    ]
    node [ 
        id 5
        c [ str "foo" ]
    ]
    edge [
        source 2 target 5
        label "edge label"
    ]
    AttOne 1
    AttOne "x"
    AttTwo 3.14
]

# second graph is ignored
graph [ node [ ] ]
