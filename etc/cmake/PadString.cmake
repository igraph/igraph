# ------------------------------------------------------------------------------
# Macro PAD_STRING
#
# This function pads a string on the left side with a specified character to
# reach the specified length. If the string length is already long enough or
# longer, the string will not be modified.
#
# PAD_STRING(OUT_VARIABLE DESIRED_LENGTH FILL_CHAR VALUE)
#
#     OUT_VARIABLE: name of the resulting variable to create
#     DESIRED_LENGTH: desired length of the generated string
#     FILL_CHAR: character to use for padding
#     VALUE: string to pad
#
# Copyright (C) 2011 by Johannes Wienke <jwienke at techfak dot uni-bielefeld dot de>
#
# This program is free software; you can redistribute it
# and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation;
# either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# ------------------------------------------------------------------------------
FUNCTION(PAD_STRING OUT_VARIABLE DESIRED_LENGTH FILL_CHAR VALUE)
    STRING(LENGTH "${VALUE}" VALUE_LENGTH)
    MATH(EXPR REQUIRED_PADS "${DESIRED_LENGTH} - ${VALUE_LENGTH}")
    SET(PAD ${VALUE})
    IF(REQUIRED_PADS GREATER 0)
        MATH(EXPR REQUIRED_MINUS_ONE "${REQUIRED_PADS} - 1")
        FOREACH(FOO RANGE ${REQUIRED_MINUS_ONE})
            SET(PAD "${FILL_CHAR}${PAD}")
        ENDFOREACH()
    ENDIF()
    SET(${OUT_VARIABLE} "${PAD}" PARENT_SCOPE)
ENDFUNCTION()

