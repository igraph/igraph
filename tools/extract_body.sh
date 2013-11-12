#!/bin/sh

sed -n '1,/^<body[ >]/!{ /<\/body>/,/^<body[ >]/!p; }'

