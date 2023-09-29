# comment lines are ignored
Version 1
graph [
    directed 1
    node [
        id 1
        a [ b 1 ]
        b +2
        c 5.6
        graphics [ x 0 y 0 ]
    ]
    node [
        a 1
        id 2
        a 2
        b "asd"
        d -Inf
        nan "foo"
    ]
    edge [
        source 5 target 1
        weight -0.45
        label -0.45
        num1 NAN
    ]
    node [
        id 5
        c [ str "bar" ]
    ]
    edge [
        source 2 target 5
        label "Tom &amp; Jerry&apos;s &quot;friendship&quot;"
# The 'sou_rce' attribute will be ignored by the GML writer as it would conflict with
# 'source' after stripping the '_' character.
        sou_rce 1.0
        num1 +inF
    ]
    AttOne 1
    AttOne "x"
    AttTwo 3.14
]

# second graph is ignored
graph [ node [ ] ]
