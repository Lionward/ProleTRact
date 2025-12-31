#!/usr/bin/env python3
"""
ProleTRact CLI - Command-line interface for starting the application
"""
import sys
import os
import subprocess
import signal
import time
import webbrowser
from pathlib import Path
from typing import Optional

# Get the package directory
# When installed, __file__ points to site-packages/proletract/cli/main.py
# When in development, it points to the source directory
_package_root = Path(__file__).parent.parent
if (_package_root / "backend").exists():
    # Development mode - backend is at proletract/backend
    BACKEND_DIR = _package_root / "backend"
    FRONTEND_DIR = _package_root / "frontend"
    PACKAGE_DIR = _package_root.parent
else:
    # Installed mode - find the package data
    import proletract
    _package_path = Path(proletract.__file__).parent
    BACKEND_DIR = _package_path / "backend"
    FRONTEND_DIR = _package_path / "frontend"
    PACKAGE_DIR = _package_path.parent


def find_node():
    """Find Node.js executable"""
    try:
        result = subprocess.run(["which", "node"], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
    except:
        pass
    
    # Try common locations
    for path in ["/usr/bin/node", "/usr/local/bin/node", "node"]:
        if os.path.exists(path) or path == "node":
            try:
                subprocess.run([path, "--version"], capture_output=True, check=True)
                return path
            except:
                continue
    return None


def is_port_in_use(port: int, host: str = "127.0.0.1") -> bool:
    """Check if a port is already in use"""
    import socket
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.settimeout(1)
            result = s.connect_ex((host, port))
            return result == 0
    except:
        return False


def get_port_process(port: int) -> Optional[str]:
    """Get the process using a port (if any)"""
    try:
        # Try lsof first (Linux/macOS)
        result = subprocess.run(
            ["lsof", "-ti", f":{port}"],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            pid = result.stdout.strip()
            if pid:
                # Get process name
                try:
                    ps_result = subprocess.run(
                        ["ps", "-p", pid, "-o", "comm="],
                        capture_output=True,
                        text=True
                    )
                    if ps_result.returncode == 0:
                        return f"{ps_result.stdout.strip()} (PID: {pid})"
                    return f"PID: {pid}"
                except:
                    return f"PID: {pid}"
    except FileNotFoundError:
        pass
    
    try:
        # Try netstat (alternative)
        result = subprocess.run(
            ["netstat", "-tuln"],
            capture_output=True,
            text=True
        )
        if f":{port}" in result.stdout:
            return "Unknown process"
    except:
        pass
    
    return None


def check_dependencies():
    """Check if all dependencies are installed"""
    issues = []
    
    # Check Python dependencies
    try:
        import fastapi
        import uvicorn
        import pysam
    except ImportError as e:
        issues.append(f"Missing Python dependency: {e.name}")
    
    # Check Node.js
    node_path = find_node()
    if not node_path:
        issues.append("Node.js not found. Please install Node.js 16+")
    
    # Check npm
    try:
        subprocess.run(["npm", "--version"], capture_output=True, check=True)
    except:
        issues.append("npm not found. Please install npm")
    
    # Check if frontend dependencies are installed
    if FRONTEND_DIR.exists():
        node_modules = FRONTEND_DIR / "node_modules"
        if not node_modules.exists():
            issues.append("Frontend dependencies not installed. Run 'npm install' in frontend directory")
    
    return issues


def install_frontend_deps():
    """Install frontend dependencies"""
    if not FRONTEND_DIR.exists():
        print("❌ Frontend directory not found!")
        return False
    
    print("📦 Installing frontend dependencies...")
    try:
        result = subprocess.run(
            ["npm", "install", "--legacy-peer-deps"],
            cwd=FRONTEND_DIR,
            check=True
        )
        print("✅ Frontend dependencies installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install frontend dependencies: {e}")
        return False
    except FileNotFoundError:
        print("❌ npm not found. Please install Node.js and npm")
        return False


def start_backend(port: int = 8502, host: str = "127.0.0.1", reload: bool = True, kill_existing: bool = False):
    """Start the FastAPI backend server"""
    if not BACKEND_DIR.exists():
        print("❌ Backend directory not found!")
        return None
    
    # Check if port is in use
    if is_port_in_use(port, host):
        process_info = get_port_process(port)
        if kill_existing:
            print(f"⚠️  Port {port} is in use by {process_info or 'unknown process'}")
            print(f"   Attempting to free port {port}...")
            try:
                subprocess.run(["lsof", "-ti", f":{port}", "|", "xargs", "kill", "-9"], 
                             shell=True, check=False)
                time.sleep(1)
            except:
                pass
        else:
            print(f"❌ Port {port} is already in use!")
            if process_info:
                print(f"   Process: {process_info}")
            print(f"   Use --backend-port to specify a different port")
            print(f"   Or kill the process: lsof -ti:{port} | xargs kill -9")
            return None
    
    # Import here to avoid issues if not installed
    try:
        import uvicorn
        from proletract.backend import main
    except ImportError as e:
        print(f"❌ Failed to import backend: {e}")
        print("   Make sure all Python dependencies are installed:")
        print("   pip install -r requirements.txt")
        return None
    
    print(f"🚀 Starting backend server on http://{host}:{port}")
    
    # Start uvicorn in a subprocess
    cmd = [
        sys.executable, "-m", "uvicorn",
        "proletract.backend.main:app",
        "--host", host,
        "--port", str(port),
    ]
    
    if reload:
        cmd.append("--reload")
    
    process = subprocess.Popen(
        cmd,
        cwd=PACKAGE_DIR,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    
    # Wait a moment to check if it started successfully
    time.sleep(2)
    if process.poll() is not None:
        stdout, stderr = process.communicate()
        error_msg = stderr.decode() if stderr else stdout.decode() if stdout else "Unknown error"
        print(f"❌ Backend failed to start:")
        # Show relevant error message
        if "Address already in use" in error_msg or "Errno 98" in error_msg:
            print(f"   Port {port} is already in use")
            process_info = get_port_process(port)
            if process_info:
                print(f"   Process using port: {process_info}")
            print(f"   Use --backend-port to specify a different port")
        else:
            print(error_msg[:500])  # Print first 500 chars
        return None
    
    return process


def start_frontend(port: int = 3000, install_deps: bool = False, use_build: bool = True, kill_existing: bool = False):
    """Start the React frontend development server or serve built app"""
    if not FRONTEND_DIR.exists():
        print("❌ Frontend directory not found!")
        return None
    
    original_port = port
    
    # Check if port is in use
    if is_port_in_use(port):
        process_info = get_port_process(port)
        if kill_existing:
            print(f"⚠️  Port {port} is in use by {process_info or 'unknown process'}")
            print(f"   Attempting to free port {port}...")
            try:
                subprocess.run(["lsof", "-ti", f":{port}", "|", "xargs", "kill", "-9"], 
                             shell=True, check=False)
                time.sleep(1)
                # Check again after killing
                if is_port_in_use(port):
                    print(f"   Port {port} still in use, trying alternative...")
                    kill_existing = False  # Fall through to find alternative
                else:
                    print(f"   Port {port} is now free")
            except:
                pass
        
        if is_port_in_use(port) and not kill_existing:
            print(f"⚠️  Port {port} is already in use!")
            if process_info:
                print(f"   Process: {process_info}")
            print(f"   Use --port to specify a different port")
            print(f"   Or kill the process: lsof -ti:{port} | xargs kill -9")
            # Try to find an available port
            for test_port in range(port + 1, port + 10):
                if not is_port_in_use(test_port):
                    print(f"   Trying alternative port {test_port}...")
                    port = test_port
                    break
            else:
                print("   No available ports found. Please free port 3000 or specify --port")
                return None
    
    # Check if built version exists
    build_dir = FRONTEND_DIR / "build"
    if use_build and build_dir.exists() and (build_dir / "index.html").exists():
        # Serve the built React app using a simple HTTP server
        print(f"🚀 Serving built frontend on http://localhost:{port}")
        try:
            import http.server
            import socketserver
            
            class Handler(http.server.SimpleHTTPRequestHandler):
                def __init__(self, *args, **kwargs):
                    super().__init__(*args, directory=str(build_dir), **kwargs)
            
            httpd = socketserver.TCPServer(("", port), Handler)
            # Run in a thread
            import threading
            thread = threading.Thread(target=httpd.serve_forever, daemon=True)
            thread.start()
            print(f"✅ Frontend served from build directory")
            return httpd  # Return the server object
        except OSError as e:
            if "Address already in use" in str(e) or "Errno 98" in str(e):
                print(f"⚠️  Port {port} is still in use, falling back to development server...")
            else:
                print(f"⚠️  Failed to serve built app: {e}")
                print("   Falling back to development server...")
            use_build = False
        except Exception as e:
            print(f"⚠️  Failed to serve built app: {e}")
            print("   Falling back to development server...")
            use_build = False
    
    # Fall back to development server
    node_path = find_node()
    if not node_path:
        print("❌ Node.js not found!")
        return None
    
    # Check if dependencies are installed
    node_modules = FRONTEND_DIR / "node_modules"
    if not node_modules.exists():
        if install_deps:
            if not install_frontend_deps():
                return None
        else:
            print("❌ Frontend dependencies not installed!")
            print("   Run with --install-deps flag to install automatically")
            return None
    
    print(f"🚀 Starting frontend development server on http://localhost:{port}")
    
    # Set environment variables for React
    env = os.environ.copy()
    env["HOST"] = "localhost"
    env["PORT"] = str(port)
    env["DANGEROUSLY_DISABLE_HOST_CHECK"] = "true"
    
    process = subprocess.Popen(
        ["npm", "start"],
        cwd=FRONTEND_DIR,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    
    return process


def main():
    """Main CLI entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="ProleTRact - Tandem Repeat Visualization Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  proletract                    # Start both servers with default settings
  proletract --backend-only     # Start only the backend server
  proletract --frontend-only    # Start only the frontend server
  proletract --port 8080        # Use custom port for frontend
  proletract --backend-port 9000 # Use custom port for backend
  proletract --install-deps     # Install frontend dependencies automatically
  proletract --no-browser        # Don't open browser automatically
        """
    )
    
    parser.add_argument(
        "--backend-only",
        action="store_true",
        help="Start only the backend server"
    )
    parser.add_argument(
        "--frontend-only",
        action="store_true",
        help="Start only the frontend server"
    )
    parser.add_argument(
        "--port",
        type=int,
        default=3000,
        help="Frontend server port (default: 3000)"
    )
    parser.add_argument(
        "--backend-port",
        type=int,
        default=8502,
        help="Backend server port (default: 8502)"
    )
    parser.add_argument(
        "--host",
        type=str,
        default="127.0.0.1",
        help="Backend server host (default: 127.0.0.1)"
    )
    parser.add_argument(
        "--install-deps",
        action="store_true",
        help="Install frontend dependencies if missing"
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Don't open browser automatically"
    )
    parser.add_argument(
        "--no-reload",
        action="store_true",
        help="Disable auto-reload for backend (production mode)"
    )
    parser.add_argument(
        "--check-deps",
        action="store_true",
        help="Check dependencies and exit"
    )
    parser.add_argument(
        "--kill-existing",
        action="store_true",
        help="Kill processes using the required ports if they're in use"
    )
    
    args = parser.parse_args()
    
    # Check dependencies
    if args.check_deps:
        print("🔍 Checking dependencies...")
        issues = check_dependencies()
        if issues:
            print("\n❌ Issues found:")
            for issue in issues:
                print(f"   - {issue}")
            sys.exit(1)
        else:
            print("✅ All dependencies are installed")
            sys.exit(0)
    
    # Check dependencies before starting
    issues = check_dependencies()
    if issues and not args.install_deps:
        print("⚠️  Some dependencies are missing:")
        for issue in issues:
            print(f"   - {issue}")
        print("\n💡 Tip: Use --install-deps to install frontend dependencies automatically")
        print("   Or run: proletract --check-deps for details")
        if "Frontend dependencies" not in str(issues):
            # Only exit if critical dependencies are missing
            pass
    
    processes = []
    
    def cleanup(signum=None, frame=None):
        """Cleanup function to stop all processes"""
        print("\n\n🛑 Shutting down servers...")
        for process in processes:
            if process:
                # Check if it's a subprocess
                if hasattr(process, 'poll'):
                    if process.poll() is None:
                        process.terminate()
                        try:
                            process.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            process.kill()
                # Check if it's an HTTP server
                elif hasattr(process, 'shutdown'):
                    try:
                        process.shutdown()
                    except:
                        pass
        sys.exit(0)
    
    # Register signal handlers
    signal.signal(signal.SIGINT, cleanup)
    signal.signal(signal.SIGTERM, cleanup)
    
    try:
        # Start backend
        if not args.frontend_only:
            backend_process = start_backend(
                port=args.backend_port,
                host=args.host,
                reload=not args.no_reload,
                kill_existing=args.kill_existing
            )
            if backend_process:
                processes.append(backend_process)
                print(f"✅ Backend running at http://{args.host}:{args.backend_port}")
                print(f"   API docs: http://{args.host}:{args.backend_port}/docs")
            else:
                if not args.backend_only:
                    print("⚠️  Backend failed to start, continuing with frontend only...")
        
        # Start frontend
        if not args.backend_only:
            time.sleep(2)  # Give backend time to start
            frontend_process = start_frontend(
                port=args.port,
                install_deps=args.install_deps,
                kill_existing=args.kill_existing
            )
            if frontend_process:
                processes.append(frontend_process)
                print(f"✅ Frontend running at http://localhost:{args.port}")
            else:
                print("⚠️  Frontend failed to start")
        
        # Open browser
        if not args.no_browser and not args.backend_only:
            time.sleep(3)  # Wait for servers to be ready
            url = f"http://localhost:{args.port}"
            print(f"\n🌐 Opening browser at {url}")
            webbrowser.open(url)
        
        print("\n" + "="*50)
        print("ProleTRact is running!")
        if not args.backend_only:
            print(f"   Frontend: http://localhost:{args.port}")
        if not args.frontend_only:
            print(f"   Backend:  http://{args.host}:{args.backend_port}")
            print(f"   API Docs: http://{args.host}:{args.backend_port}/docs")
        print("="*50)
        print("\nPress Ctrl+C to stop\n")
        
        # Wait for processes
        for process in processes:
            if process:
                process.wait()
    
    except KeyboardInterrupt:
        cleanup()
    except Exception as e:
        print(f"\n❌ Error: {e}")
        cleanup()
        sys.exit(1)


if __name__ == "__main__":
    main()

