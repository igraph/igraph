
#   IGraph R package
#   Copyright (C) 2009  Gabor Csardi <csardi@rmki.kfki.hu>
#   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

proc render {win href} {
    global tkigraph_help_root
    set Url $tkigraph_help_root/$href
    $win configure -state normal
    HMreset_win $win
    HMparse_html [get_html $Url] "HMrender $win"
    $win configure -state disabled
    update
}

proc HMlink_callback {win href} {
    render $win $href
}

proc get_html {file} {
    global tkigraph_help_root
    if {[catch {set fd [open $file]} msg]} {
	return "
                        <title>Bad file $file</title>
                        <h1>Error reading $file</h1><p>
                        $msg<hr>
                        <a href=$tkigraph_help_root>Go home</a>
                "
    }
    set result [read $fd]
    close $fd
    return $result
}
