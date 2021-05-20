# delete the header produced by f2c
/\/\*  -- trans.*/,/\*\//{
d
}

# extract the first line of including f2c.h
/#include/,/^$/{
w temp/header1
d
}

# possible local constants produced by f2c
/\/\* Table of constant values \*\//{
s/^/    /
w temp/header2
d
}
/^static.*=.*/,/^$/{
s/^/    /
w temp/header2
d
}

# matches /* Subroutine */..._(  or /* Complex */..._(
/^\/\* .*_(/,/^\{/{
w temp/header3
/^\{/!{
d
}
}
# matches any function declaration line
/^[a-zA-Z].*_(/,/^\{/{
w temp/header3
/^\{/!{
d
}
}

/^\{/,/\/\*.*LAPACK/{
/\/\*.*LAPACK/!{
/^$/d
/^\{/d
w temp/prologue
d
}
}

/\/\*.*LAPACK/,/\*\//{
w temp/comment
d
}


w temp/code
