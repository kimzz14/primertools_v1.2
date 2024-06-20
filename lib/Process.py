import subprocess, threading

class PROC:
    def __init__(self, command, pid = -1, status = 'ready'):
        self.process = None
        self.command = command
        self.pid = pid
        self.status = status
        self.thread = None
        self.stdout = []
        self.stderr = []
    
    def run(self):
        global devnull
        self.process = subprocess.Popen(self.command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        self.pid = self.process.pid
        self.status = 'running'

        self.thread = threading.Thread(target=self.read_process, args=())
        self.thread.start()
        return self

    def read_process(self):
        self.stdout = self.process.stdout.read().decode('utf-8')
        self.stderr = self.process.stderr.read().decode('utf-8')
    def wait(self):
        self.process.wait()
        self.thread.join()
        return self

    def get_status(self):
        if self.status == 'complete':
            pass
        elif self.status == 'ready':
            pass
        elif self.status == 'running':
            ps_process = PROC('ps --pid ' +  str(self.pid) + ' -o command')
            ps_process.run().wait()
            ps_command = ps_process.stdout.split('\n')[1]
            if ps_command == '':
                self.status = 'complete'
            elif ps_command.find('defunct') != -1:
                self.status = 'complete'
        else:
            print('bug!')
        return self.status


if __name__ == "__main__":
    process = PROC('python test.py')
    process.run()
    print(process.pid)
    print(process.get_status())

    print(process.stdout)
    print(process.stderr)