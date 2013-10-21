#! /bin/sh

git describe HEAD --tags | rev | sed 's/g-/./' | sed 's/-/+/' | rev
