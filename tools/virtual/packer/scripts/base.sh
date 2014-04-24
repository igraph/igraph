#!/bin/sh -x

#Turn off and disable ufw
echo "Stopping ufw..."
service ufw stop
echo "Disabling ufw..."
ufw disable