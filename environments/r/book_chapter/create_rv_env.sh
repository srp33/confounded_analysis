#!/bin/bash
#
# Script to create and update a self-contained R project using rig and rv.

# --- Configuration ---
R_VERSION="4.4"
PROJECT_NAME="book_chapter"

echo "--- Checking if R version $R_VERSION has been installed by rig ---"
if ! rig ls | grep -q "^$R_VERSION"; then
  echo "R $R_VERSION not found. Installing with rig..."
  rig add $R_VERSION
else
  echo "R $R_VERSION is already installed."
fi

# --- Sync the Environment ---
echo "\n--- Planning environment (dry run) ---"
rv plan

echo "\n--- Syncing environment (installing packages) ---"
rv sync

echo "\n✅ Done. Your new project is ready in ./$PROJECT_NAME"