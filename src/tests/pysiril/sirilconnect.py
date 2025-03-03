import os
import sys
import threading
import time
import queue
import platform
import re
import subprocess

class SirilPipeController:
    """
    A cross-platform controller class for interacting with Siril CLI via named pipes.
    Works on both Windows and Linux/macOS platforms with improved command completion detection.
    """
    
    def __init__(self, pipe_in_name="siril_command.in", pipe_out_name="siril_command.out"):
        self.pipe_in_name = pipe_in_name
        self.pipe_out_name = pipe_out_name
        self.in_pipe = None
        self.out_pipe = None
        self.output_queue = queue.Queue()
        self.output_thread = None
        self.running = False
        self.last_command_output = []
        self.command_finished = threading.Event()
        self.is_windows = platform.system() == "Windows"
        
        # Pattern matching for command prompt
        # # Adjust this regex based on Siril's actual prompt pattern
        # self.prompt_pattern = re.compile(r'(siril|SiriL)(:|\>|\$)')
        self.stable_output_timeout = 10.  # Time with no new output to consider command complete
        self.last_output_time = 0
        
        # Siril process
        self.siril_process = None
    
    def connect(self):
        """Connect to both input and output pipes."""
        try:
            if self.is_windows:
                return self._connect_windows()
            else:
                return self._connect_unix()
        except Exception as e:
            print(f"Error connecting to pipes: {e}")
            self.disconnect()
            return False
    
    def _connect_windows(self):
        """Windows-specific pipe connection logic using the Win32 API."""
        import win32pipe
        import win32file
        
        # Connect to the input pipe (for sending commands)
        self.in_pipe = win32file.CreateFile(
            r'\\.\pipe\{}'.format(self.pipe_in_name),
            win32file.GENERIC_WRITE,
            0,
            None,
            win32file.OPEN_EXISTING,
            0,
            None
        )
        
        # Connect to the output pipe (for receiving responses)
        self.out_pipe = win32file.CreateFile(
            r'\\.\pipe\{}'.format(self.pipe_out_name),
            win32file.GENERIC_READ,
            0,
            None,
            win32file.OPEN_EXISTING,
            0,
            None
        )
        
        # Start the output reader thread
        self.running = True
        self.output_thread = threading.Thread(target=self._read_output_windows)
        self.output_thread.daemon = True
        self.output_thread.start()
        
        return True
    
    def _connect_unix(self):
        """Unix-specific pipe connection logic using file-based pipes."""
        # Get path to the pipes - typically in /tmp on Linux
        pipe_in_path = f"/tmp/{self.pipe_in_name}"
        pipe_out_path = f"/tmp/{self.pipe_out_name}"
        
        # Check if pipes exist
        if not os.path.exists(pipe_in_path) or not os.path.exists(pipe_out_path):
            raise FileNotFoundError(f"Pipes not found at {pipe_in_path} and/or {pipe_out_path}")
        
        # Open the pipes
        self.in_pipe = open(pipe_in_path, 'w')
        self.out_pipe = open(pipe_out_path, 'r')
        
        # Start the output reader thread
        self.running = True
        self.output_thread = threading.Thread(target=self._read_output_unix)
        self.output_thread.daemon = True
        self.output_thread.start()
        
        return True
    
    def disconnect(self):
        """Disconnect from both pipes and clean up resources."""
        self.running = False
        
        if self.output_thread and self.output_thread.is_alive():
            self.output_thread.join(timeout=1.0)

        # Terminate Siril process if it was started by us
        if self.siril_process:
            try:
                self.siril_process.terminate()
                self.siril_process.wait(timeout=5)
            except:
                if self.is_windows:
                    os.system(f"taskkill /F /PID {self.siril_process.pid}")
                else:
                    os.system(f"kill -9 {self.siril_process.pid}")
            self.siril_process = None
    
    def _check_command_completion(self, line):
        """
        Check if the line indicates command completion based on multiple strategies.
        Returns True if the command is likely complete.
        """
        
        if 'status: success' in line or 'status: error' in line:
            return True
            
        # Update the last output time
        self.last_output_time = time.time()
        return False
    
    def _read_output_windows(self):
        """Windows-specific thread function to read from output pipe."""
        import win32file
        
        buffer = ""
        self.last_output_time = time.time()
        
        while self.running:
            try:
                # Check for stable output timeout (no new output for a while)
                if (not self.command_finished.is_set() and 
                    self.last_command_output and 
                    time.time() - self.last_output_time > self.stable_output_timeout):
                    self.command_finished.set()
                
                # Non-blocking read with timeout
                hr, data = win32file.ReadFile(self.out_pipe, 4096)
                if data:
                    text = data.decode('utf-8', errors='replace')
                    buffer += text
                    
                    # Process complete lines
                    lines = buffer.split('\n')
                    buffer = lines.pop()  # Keep the last incomplete line
                    
                    for line in lines:
                        line = line.rstrip('\r')
                        if line:
                            self.output_queue.put(line)
                            self.last_command_output.append(line)
                            
                            # Check if line indicates command completion
                            if self._check_command_completion(line):
                                self.command_finished.set()
            except Exception as e:
                if not isinstance(e, win32file.error):  # Ignore timeout errors
                    print(f"Error reading from output pipe: {e}")
                    break
            
            # Small sleep to prevent CPU hogging
            time.sleep(0.01)
    
    def _read_output_unix(self):
        """Unix-specific thread function to read from output pipe."""
        self.last_output_time = time.time()
        
        # Make the pipe non-blocking
        import fcntl
        import os
        fd = self.out_pipe.fileno()
        fl = fcntl.fcntl(fd, fcntl.F_GETFL)
        fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
        
        buffer = ""
        while self.running:
            try:
                # Check for stable output timeout (no new output for a while)
                if (not self.command_finished.is_set() and 
                    self.last_command_output and 
                    time.time() - self.last_output_time > self.stable_output_timeout):
                    self.command_finished.set()
                
                # Try to read a chunk of data
                chunk = self.out_pipe.read(4096)
                if chunk:
                    buffer += chunk
                    
                    # Process complete lines
                    lines = buffer.split('\n')
                    buffer = lines.pop()  # Keep the last incomplete line
                    
                    for line in lines:
                        line = line.rstrip('\r')
                        if line:
                            self.output_queue.put(line)
                            self.last_command_output.append(line)
                            
                            # Check if line indicates command completion
                            if self._check_command_completion(line):
                                self.command_finished.set()
            except (IOError, BlockingIOError):
                # Expected for non-blocking IO when no data is available
                pass
            except Exception as e:
                print(f"Error reading from output pipe: {e}")
                break
            
            # Small sleep to prevent CPU hogging
            time.sleep(0.01)
    
    def send_command(self, command, timeout=10.0, debug=False):
        """
        Send a single command and wait for its completion.
        
        Args:
            command: The command to send
            timeout: Maximum time to wait for a response in seconds
            debug: Whether to print debug information
            
        Returns:
            A list of output lines received in response to the command
        """
        if not self.in_pipe or not self.out_pipe:
            raise RuntimeError("Not connected to pipes. Call connect() first.")
        
        if debug:
            print(f"Sending command: {command}")
        
        # Clear previous command output and reset event
        self.last_command_output = []
        self.command_finished.clear()
        self.last_output_time = time.time()
        
        # Send the command based on platform
        if self.is_windows:
            self._send_command_windows(command)
        else:
            self._send_command_unix(command)
        
        # Wait for the command to complete
        start_time = time.time()
        command_completed = False
        
        while not command_completed and (time.time() - start_time < timeout):
            # Check if command has completed
            if self.command_finished.is_set():
                command_completed = True
                break
                
            # Wait a bit before checking again
            self.command_finished.wait(0.1)
            
            # Also consider command complete if no output for a while
            if self.last_command_output and (time.time() - self.last_output_time > self.stable_output_timeout):
                command_completed = True
                if debug:
                    print("Command considered complete due to stable output")
                break
                
        if not command_completed:
            if debug:
                print(f"Warning: Command timed out after {timeout} seconds: {command}")
        elif debug:
            print(f"Command completed in {time.time() - start_time:.2f} seconds")
            
        return self.last_command_output.copy()
    
    def _send_command_windows(self, command):
        """Windows-specific logic to send a command."""
        import win32file
        command_bytes = (command + "\n").encode('utf-8')
        win32file.WriteFile(self.in_pipe, command_bytes)
    
    def _send_command_unix(self, command):
        """Unix-specific logic to send a command."""
        self.in_pipe.write(command + "\n")
        self.in_pipe.flush()
    
    def send_commands(self, commands, timeout_per_command=10.0, debug=False):
        """
        Send multiple commands in sequence and collect all outputs.
        
        Args:
            commands: List of commands to send
            timeout_per_command: Maximum time to wait for each command
            debug: Whether to print debug information
            
        Returns:
            A dictionary mapping commands to their respective outputs
        """
        results = {}
        
        for cmd in commands:
            output = self.send_command(cmd, timeout=timeout_per_command, debug=debug)
            results[cmd] = output
            
        return results
    
    def get_pending_output(self, max_lines=None):
        """
        Get any pending output from the queue without waiting.
        
        Args:
            max_lines: Maximum number of lines to retrieve (None for all available)
            
        Returns:
            List of output lines
        """
        lines = []
        try:
            while max_lines is None or len(lines) < max_lines:
                lines.append(self.output_queue.get_nowait())
                self.output_queue.task_done()
        except queue.Empty:
            pass  # No more items
            
        return lines

    def launch_siril(self, executable_path=None, arguments=None):
        """
        Launch Siril CLI application.
        
        Args:
            executable_path: Path to Siril executable
            arguments: List of command-line arguments
            
        Returns:
            True if launched successfully, False otherwise
        """
        # Default paths based on platform
        if executable_path is None:
            if self.is_windows:
                executable_path = r"C:\msys64\mingw64\bin\siril-cli.exe"
            else:
                executable_path = "siril-cli"
                
        # Default arguments
        if arguments is None:
            arguments = ["-p", "-o"]  # Use pipes for communication and offline mode
            
        # Build the command
        cmd = [executable_path] + arguments
        
        try:
            # Launch the process
            print(f"Launching Siril: {' '.join(cmd)}")
            self.siril_process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                creationflags=subprocess.CREATE_NO_WINDOW if self.is_windows else 0
            )
            
            # Give it a moment to start up and create the pipes
            time.sleep(1)
            
            return True
        except Exception as e:
            print(f"Error launching Siril: {e}")
            return False


# Main function
def main():
    """
    Main function that demonstrates launching Siril and interacting with it.
    """
    # Import required platform-specific modules for Windows
    if platform.system() == "Windows":
        try:
            import win32pipe
            import win32file
        except ImportError:
            print("On Windows, this module requires the pywin32 package.")
            print("Please install it with: pip install pywin32")
            sys.exit(1)
        sirilexe = r"C:\msys64\mingw64\bin\siril-cli.exe"
    else:
        sirilexe = "siril-cli"
    
    # Create the controller
    controller = SirilPipeController()
    
    # Launch Siril
    if not controller.launch_siril(sirilexe, ["-p"]):
        print("Failed to launch Siril.")
        return
    
    # Wait for pipes to be created
    print("Waiting for Siril to initialize...")
    time.sleep(2)
    
    # Connect to the pipes
    if not controller.connect():
        print("Failed to connect to Siril pipes.")
        return
    
    try:
        print("Connected to Siril successfully!")
        
        # Enable debug mode
        debug_mode = True
        
        # Example: Run some basic commands
        print("\nChecking Siril version...")
        output = controller.send_command("requires 1.3.6", debug=debug_mode)
        print(f"Received {len(output)} lines of output:")
        for line in output:
            print(f"  > {line}")
        
        # Example: Run multiple commands
        print("\nRunning multiple commands...")
        path = os.path.dirname(__file__)
        sigma = 0.7
        roundness = 0.75
        convergence = 3
        commands = [
            f"cd {path}",
            "setext fit",
            "load star_test1",
            "findstar -out=starlist1.lst",
            f"setfindstar -sigma={sigma} -roundness={roundness} -convergence={convergence}",
            "findstar -out=starlist2.lst",
            "close",
        ]
        
        results = controller.send_commands(commands, debug=debug_mode)
        
        for cmd, out in results.items():
            print(f"\nCommand: {cmd}")
            print(f"Received {len(out)} lines:")
            for line in out:
                print(f"  > {line}")

    finally:
        # Clean up resources
        controller.disconnect()
        print("Disconnected from Siril and terminated process.")


# Run the main function if script is executed directly
if __name__ == "__main__":
    main()