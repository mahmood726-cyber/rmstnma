#!/bin/bash
# powerNMA Startup Script

set -e

echo "Starting powerNMA services..."

# Start Ollama service in background
if command -v ollama &> /dev/null; then
    echo "Starting Ollama service..."
    ollama serve &
    sleep 5

    # Pull default model
    if [ ! -d "/root/.ollama/models" ]; then
        echo "Pulling llama2 model..."
        ollama pull llama2
    fi
fi

# Start REST API in background (if enabled)
if [ "$POWERNMA_ENABLE_API" = "true" ]; then
    echo "Starting REST API..."
    R -e "library(powerNMA); launch_powernma_api(host='0.0.0.0', port=8000)" &
fi

# Start Shiny Server
echo "Starting Shiny Server..."
exec shiny-server
