#! /bin/bash
# download data from LNE-OP repository (consortium-only access)
# which  needs to be setup un in advance using rclone as obs-nextcloud2025
rclone sync obs-nextcloud2025: ./Data --progress

