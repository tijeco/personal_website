# Use an official TeX Live image as a base
FROM texlive/texlive:latest

# Install latexmk
RUN apt-get update && apt-get install -y \
    latexmk \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /workspace

# Set the default command to run latexmk