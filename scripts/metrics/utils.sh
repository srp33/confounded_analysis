archive_file() {
  local file_path="$1"
  if [ -f "${file_path}" ]; then
    local mod_date
    mod_date=$(date -r "${file_path}" +%Y-%m-%d-%h)
    local filename
    filename=$(basename -- "${file_path}")
    cp "${file_path}" "outputs/metrics/archive/${mod_date}_${filename}"
  fi
}