#!/usr/bin/env python3
import subprocess
import fileinput
import tarfile
import shutil
import os
import sys

#from https://stackoverflow.com/questions/14697629/running-a-bat-file-though-python-in-current-process
def init_vsvars():
    vswhere_path = r"%ProgramFiles(x86)%/Microsoft Visual Studio/Installer/vswhere.exe"
    vswhere_path = os.path.expandvars(vswhere_path)
    if not os.path.exists(vswhere_path):
        raise EnvironmentError("vswhere.exe not found at: %s", vswhere_path)

    vs_path = os.popen('"{}" -latest -property installationPath'.format(vswhere_path)).read().rstrip()
    vsvars_path = os.path.join(vs_path, "VC\\Auxiliary\\Build\\vcvars64.bat")

    output = os.popen('"{}" && set'.format(vsvars_path)).read()

    for line in output.splitlines():
        pair = line.split("=", 1)
        if(len(pair) >= 2):
            os.environ[pair[0]] = pair[1]


def get_ver():
    with open("src/common/defs.h") as f:
        for line in f.readlines():
            line = line.strip()
            if "BKC_VERSION" in line:
                return line.split("BKC_VERSION")[-1].split("\"")[-2].strip()
    print("Error: cannot read BKC version")
    sys.exit(1)

def get_os():
    if os.name == 'nt':
        return 'windows'
    elif os.name == 'posix':
        if os.uname()[0] == 'Linux':
            return 'linux'
        elif os.uname()[0] == 'Darwin':
            return 'mac'
        else:
            print("Error: unknown os", os.uname()[0])
            sys.exit(1)
    else:
        print("Error: unknown os.name", os.name)
        sys.exit(1)

def get_hardware():
    if os.name == 'nt':
        return 'x64' # TODO: do a real check and support ARM also...
    elif os.name == 'posix':
        if os.uname()[4] == 'x86_64':
            return 'x64'
        elif os.uname()[4] == 'aarch64' or os.uname()[4] == 'arm64':
            return 'arm64'
        else:
            print("Error: unknown hardware", os.uname()[4])
            sys.exit(1)
    else:
        print("Error: unknown os.name", os.name)
        sys.exit(1)

def run_cmd(cmd):    
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()

def run_cmd_get_stdout(cmd):
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return p.stdout.decode('utf-8')


if __name__ == "__main__":
    system = get_os()
    hardware = get_hardware()

    ver = get_ver()

    print(f"building\n\tVersion: {ver}\n\tOperating system: {system}\n\tHardware: {hardware}")

    run_cmd("git submodule update --init --recursive --jobs=8")
    
    make_command = "make"
    if system == "mac":
        make_command = "gmake"

    platform = "axv"

    if hardware == "arm64":
        if system == "mac":
            platform = "m1"
        else:
            platform = "arm8"

    # In general use the default g++, but not on mac where the default is just clang++
    # which is currently not supported
    cxx = "g++"
    cc = "gcc"

    if system == "mac":
        for version in [13, 12, 11, 10]:
            if shutil.which(f"g++-{version}") and shutil.which(f"gcc-{version}"):
                out = run_cmd_get_stdout(f"g++-{version} --version")
                if "gcc" in out.lower():
                    cxx = f"g++-{version}"
                else:
                    continue

                # lets check if the same version works for CC and is GNU
                out = run_cmd_get_stdout(f"gcc-{version} --version")
                if "gcc" in out.lower():
                    cc = f"gcc-{version}"
                    break
    
    if system != 'windows':
        if not "g++" in run_cmd_get_stdout(f"{cxx} --version").lower() or not "gcc" in run_cmd_get_stdout(f"{cc} --version"):
            print(f"The selected C++ compiler ({cxx}) or C compiler ({cc}) is not GNU g++/gcc.\n"
                "If you are using macOS, you may install it with Homebrew (https://brew.sh/)")
            sys.exit(1)

    if system == 'windows':
        init_vsvars()
        #run_cmd("MSBuild.exe bkc.sln /property:Configuration=Release /property:Platform=x64")
        run_cmd("devenv bkc.sln /Build \"Release|x64\"")
        #run_cmd("devenv bkc.sln /Rebuild \"Release|x64\"")

        with tarfile.open(f"bkc-{ver}.{system}.{hardware}.tar.gz", "w:gz") as tar:
            #tar.add(source_dir, arcname=os.path.basename(source_dir))
            tar.add("x64/Release/bkc.exe", arcname="bkc.exe")
            tar.add("x64/Release/bkc_dump.exe", arcname="bkc_dump.exe")
    else:
        run_cmd(f"{make_command} clean")
        run_cmd(f"{make_command} CXX={cxx} CC={cc} PLATFORM=generic STATIC_LINK=true -j")
        run_cmd(f"cd bin; tar -c * | pigz > ../bkc-{ver}.{system}.{hardware}.tar.gz; cd ..;")
        run_cmd(f"{make_command} clean")
