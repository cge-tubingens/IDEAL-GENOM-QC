import shutil
import subprocess

class ToolNotFoundError(Exception):
    pass

def is_tool_available(tool_name: str) -> bool:
    """Check whether `tool_name` is on PATH and marked as executable."""
    return shutil.which(tool_name) is not None

def get_tool_version(tool_name: str, version_arg="--version") -> str:
    """Return the version string of the external tool, if available."""
    if not is_tool_available(tool_name):
        raise ToolNotFoundError(f"Tool '{tool_name}' is not installed or not on PATH.")

    try:
        result = subprocess.run(
            [tool_name, version_arg],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=5
        )
        return result.stdout.strip()
    except Exception as e:
        return f"Failed to get version for {tool_name}: {e}"

def check_required_tools(tools: list) -> None:
    """Check a list of required tools and raise error if any are missing."""
    missing = [tool for tool in tools if not is_tool_available(tool)]
    if missing:
        raise ToolNotFoundError(f"The following required tools are missing: {', '.join(missing)}")
