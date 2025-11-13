


r_modules=(
    r/4.4.0-ncfmhh4
    r/4.5.0-xcvdvru
    r/4.5.1-gg7txi7
    r/4.5.1-hbue2wm
    r/4.5.1-5sqddv2
    r/4.5.1-264p7tz
)

for r_module in "${r_modules[@]}"; do
    (
        # module reset &> /dev/null
        module_output=$(module load "$r_module" 2>&1)
        module load "$r_module" &> /dev/null
        if [[ $module_output == *"gcc-runtime"* ]]; then
            echo "$r_module had gcc-runtime warning"
        else 
            echo "$r_module passed gcc-runtime"
        fi

        if command -v R &> /dev/null; then
            echo "$r_module made it into PATH"
        else
            echo "Failed to load R module: $r_module"
        fi

        echo ""
    )
done