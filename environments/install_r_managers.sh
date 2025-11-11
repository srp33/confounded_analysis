# Add HOME/.local/bin to PATH if not there already
case ":$PATH:" in
  *":$HOME/.local/bin:"*)
    ;;
  *)
    export PATH="$HOME/.local/bin:$PATH"
    ;;
esac

# Install r-lib/rig (The R *language* manager)
echo "--- Checking rig installation ---"
if ! command -v rig &> /dev/null; then
    echo "Installing rig..."
    cd ~
    
    # The rig tarball is designed to be extracted directly into ~/.local
    # This will place 'bin/rig' into '~/.local/bin/rig'
    # and 'share/...' into '~/.local/share/'
    curl -Ls https://github.com/r-lib/rig/releases/download/latest/rig-linux-$(arch)-latest.tar.gz | \
      tar xz -C $HOME/.local

else
    echo "rig is already installed."
fi

# For bash, rig wants access to bash-completion
# This checks if the file we are about to install already exists
if [ ! -f "$HOME/.local/share/bash-completion/bash_completion" ]; then
    echo "--- bash-completion not found. Installing... ---"
    # Create a directory for source code
    mkdir -p ~/.local/src
    cd ~/.local/src
    curl -Ls https://github.com/scop/bash-completion/releases/download/2.17.0/bash-completion-2.17.0.tar.xz -o bash-completion.tar.xz
    tar -xf bash-completion.tar.xz
    cd bash-completion-2.17.0
    ./configure --prefix=$HOME/.local
    make
    make install

    # Add the activation line to ~/.bashrc (Check if the line already exists first)
    COMPLETION_LINE="source $HOME/.local/share/bash-completion/bash_completion"
    if ! grep -Fxq "$COMPLETION_LINE" ~/.bashrc; then
        echo "" >> ~/.bashrc
        echo "# Load local bash-completion" >> ~/.bashrc
        echo "$COMPLETION_LINE" >> ~/.bashrc
    fi
    source $HOME/.local/share/bash-completion/bash_completion
else
    echo "--- bash-completion is already installed. Skipping. ---"
fi



cd ~
# Install A2-ai/rv (The R *package* manager) 
if ! command -v rv &> /dev/null; then
    echo "Installing rv..."
    # Uses a smart installation script
    curl -sSL https://raw.githubusercontent.com/A2-ai/rv/refs/heads/main/scripts/install.sh | bash
else
    echo "rv is already installed."
fi


echo "--- Verifying installations ---"
rv --version
rig --version