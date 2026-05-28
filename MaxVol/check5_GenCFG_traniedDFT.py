import os
import subprocess
import sys
import time
import shutil
import re
import glob
import tarfile

"""
trained_dft.cfg generation from All training steps
"""

# 환경 변수 설정
mlip_path  = "/home/sulee/MLIP/bin"
local_CCpy = "/home/sulee/CCpy/CCpy/VASP"
root = os.path.basename(os.getcwd())
file_path = "./MLIP5.py"
elist = subprocess.getoutput(f"grep elist {file_path} | head -1 | awk '{{print $3, $4}}'").split()
#print(elist)

# 가상 환경 확인
try:
    env_chk = subprocess.run(
        ["which", "CCpyVASPAnal.py"], capture_output=True, text=True, check=True
    ).stdout.count("conda")
    if env_chk != 1:
        print("Please change to Virtual ENV for MLIP!!!")
        sys.exit()
except subprocess.CalledProcessError:
    print("Error checking virtual environment. Make sure `which CCpyVASPAnal.py` works.")
    sys.exit()


# 진행률 표시 함수
def print_progress_bar(iteration, total, file_name, length=40):
    percent = f"{100 * (iteration / float(total)):.1f}"
    filled_length = int(length * iteration // total)
    bar = "█" * filled_length + "-" * (length - filled_length)
    sys.stdout.write(f"\rProcessing {file_name:30} |{bar}| {percent}% Complete  ")
    sys.stdout.flush()

def process_cfg_files(original_file, new_file, output_file):
    # CFG 파일 읽기 함수
    def read_cfg_file(file_path):
        with open(file_path, 'r') as file:
            return file.readlines()

    # 기존 CFG 데이터 파일 읽기
    original_lines = read_cfg_file(original_file)
    original_line_count = len(original_lines)

    # 새로운 CFG 데이터 파일 읽기
    new_cfg_lines = read_cfg_file(new_file)

    # 정규 표현식을 사용하여 AtomData 섹션을 찾고 두 번째 열의 값을 1로 변경
    updated_cfg_lines = []
    inside_atomdata = False
    current_line_number = 0

    for line in new_cfg_lines:
        current_line_number += 1
        if line.startswith(" AtomData"):  # AtomData 섹션 시작 확인
            inside_atomdata = True
            updated_cfg_lines.append(line)
            continue

        if inside_atomdata:
            if line.strip() == "Energy":  # AtomData 섹션 종료 확인
                inside_atomdata = False
                updated_cfg_lines.append(line)
                continue

            if current_line_number > original_line_count:  # 기존 라인 수 이후의 새로운 데이터만 수정
                parts = line.split()
                if len(parts) > 1:
                    parts[1] = "1"  # 두 번째 열을 1로 변경
                formatted_line = f"{int(parts[0]):>14} {int(parts[1]):>4} {float(parts[2]):>14.6f} {float(parts[3]):>13.6f} {float(parts[4]):>13.6f} {float(parts[5]):>12.6f} {float(parts[6]):>11.6f} {float(parts[7]):>11.6f}\n"
                updated_cfg_lines.append(formatted_line)
            else:
                updated_cfg_lines.append(line)
        else:
            updated_cfg_lines.append(line)

    # 결과를 파일로 저장
    with open(output_file, 'w') as file:
        file.writelines(updated_cfg_lines)

#    print(f"Updated CFG file saved to {output_file}")


# Function to extract con values from in.cfg
def extract_con_values(cfg_path):
    try:
        result = subprocess.run(["grep", "con", cfg_path], capture_output=True, text=True, check=True)
        return [line.split()[2] for line in result.stdout.splitlines()]  # Extract third column values from each line
    except (IndexError, subprocess.CalledProcessError):
        return []

# Function to update feature line
def update_feature_lines(cfg_file, root, con_values, d_list, con_index):
    pattern = re.compile(r"^ Feature\s+EFS_by\s+VASP")
    
    with open(cfg_file, "r") as file:
        lines = file.readlines()
    
    con_values = con_values[:len(d_list)]  # Ensure enough values are available

    updated_lines = []
    
    for line in lines:
        if pattern.match(line) and con_index < len(con_values):
            new_feature_line = f" Feature   Structure\t{con_values[con_index]}\n"
            updated_lines.append(new_feature_line)
        else:
            updated_lines.append(line)
    
    with open(cfg_file, "w") as file:
        file.writelines(updated_lines)

#    print(f"pdated Feature line in {cfg_file} → \"Feature\tStructure\t{root}-{d}\"")

def merge_out(outcar_path):
    """ 특정 OUTCAR 파일을 기준으로 error*.gz 파일을 병합하고 정리 """

    # 1. outcar_path가 있는 디렉토리에서 error*.gz 파일 찾기
    outcar_dir = os.path.dirname(outcar_path)  # OUTCAR이 있는 디렉토리
    gz_files = glob.glob(os.path.join(outcar_dir, "error*.gz"))

    # 2. error*.gz 파일이 없으면 실행하지 않음
    if not gz_files:
        os.system(f"Skipping {outcar_dir}: No 'error*.gz' files found.")
        return  # 단순히 스킵하고 종료

    for gz_file in gz_files:
        # 3. 파일 이름을 디렉토리 이름으로 사용
        base_name = os.path.splitext(os.path.basename(gz_file))[0]  # 확장자 제거
        dir_name = os.path.join(outcar_dir, base_name)  # 해당 디렉토리에 생성

        # 디렉토리 생성
        os.makedirs(dir_name, exist_ok=True)

        # 4. .gz 파일을 디렉토리로 복사
        dest_path = os.path.join(dir_name, os.path.basename(gz_file))
        shutil.copy(gz_file, dest_path)

        # 5. 압축 해제 (tar -zxf 실행)
        with tarfile.open(dest_path, "r:gz") as tar:
            tar.extractall(path=dir_name)

    # 6. 현재 디렉토리의 OUTCAR 백업
    outcar_orig_path = outcar_path + "-orig"
    if os.path.exists(outcar_path):
        shutil.copy(outcar_path, outcar_orig_path)

    # 7. cat OUTCAR-orig 생성된 directory/OUTCAR > OUTCAR
    latest_dir = max(glob.glob(os.path.join(outcar_dir, "error*")), key=os.path.getmtime)  # 최신 생성된 디렉토리 찾기
    new_outcar_path = os.path.join(latest_dir, "OUTCAR")

    if os.path.exists(new_outcar_path):
        with open(outcar_path, "w") as outcar:
            # OUTCAR-orig 파일 쓰기
            with open(outcar_orig_path, "r") as orig:
                outcar.write(orig.read())
            # 새로 생성된 OUTCAR 파일 이어쓰기
            with open(new_outcar_path, "r") as new:
                outcar.write(new.read())

    # 8. OUTCAR에서 특정 구간 삭제 (General timing ~ Voluntary context switches)
    with open(outcar_path, "r") as f:
        lines = f.readlines()

    start_keyword = "General timing and accounting informations for this job:"
    end_keyword = "Voluntary context switches:"

    filtered_lines = []
    skip = False  # 삭제 여부를 판단하는 플래그

    for line in lines:
        if start_keyword in line:
            skip = True  # 특정 구간 발견 시 삭제 시작
        if not skip:
            filtered_lines.append(line)
        if skip and end_keyword in line:  # 삭제 종료 조건
            skip = False

    # 변경된 내용을 다시 OUTCAR 파일에 저장
    with open(outcar_path, "w") as f:
        f.writelines(filtered_lines)

    if os.path.exists(dir_name):
        shutil.rmtree(dir_name)
#    os.system(f'echo " Merged OUTCAR: {outcar_path}"')

def restore_outcar(outcar_path):
    """ OUTCAR-orig가 존재하면 OUTCAR로 복원 """
    outcar_orig_path = outcar_path + "-orig"
    if os.path.exists(outcar_orig_path):
        shutil.move(outcar_orig_path, outcar_path)
#        os.system(f'echo " Restored {outcar_path} from {outcar_orig_path}"')
    else:
        os.system(f'echo " Skipping restoration: No backup found for {outcar_path}"')


# Step D-DFT
def Step_D():
    base_dir = os.getcwd()

    dirs_step = [d for d in os.listdir(base_dir) if os.path.isdir(d) and d[0].isdigit() and "00" not in d]
    tot_dirs = len(dirs_step)

    for idx, d_step in enumerate(dirs_step, start=1):
        d_dft_path = os.path.join(base_dir, d_step, "D-DFT")
        if os.path.exists(d_dft_path):
            os.chdir(d_dft_path)
        else:
            print("\n")
            shutil.move("out.cfg", "trained_dft.cfg")
            print(f"Renamed out.cfg to trained_dft.cfg (overwritten if existed).")
            # os.rename("out.cfg", "trained_dft.cfg")
            # print(f"Renamed out.cfg to trained_dft.cfg.")
            sys.exit()

        # Extract con values from in.cfg
        in_cfg_path = os.path.join(d_dft_path, "in.cfg")
        con_values = extract_con_values(in_cfg_path)
        con_index = 0

        if idx != 1:
            subprocess.run(f"cp {base_dir}/out.cfg ./", shell=True, check=False)            

        ########## OUTCAR Merge ######################
        # 9. 현재 디렉토리에서 error* 디렉토리 찾기 (error*.gz 파일이 있는 디렉토리만 선택)
        dirs_err = [
            dir for dir in sorted(glob.glob("dir_POSCAR*"), key=os.path.getmtime)
            if glob.glob(os.path.join(dir, "error*.gz"))  # 해당 디렉토리에 error*.gz 파일이 있는지 확인
        ]

        # 10. 각 필터링된 디렉토리에 대해 merge_out 실행 후 OUTCAR 복원
        for dir in dirs_err:
            gz_file = os.path.join(dir, "OUTCAR.gz")
            outcar_file = os.path.join(dir, "OUTCAR")
            
            # 먼저, dir/OUTCAR.gz 파일이 존재하면 gunzip을 이용해 압축 해제합니다.
            if os.path.exists(gz_file):
                os.system(f"gunzip {gz_file}")
            
            # 압축 해제 후, OUTCAR 파일이 존재하면 병합 작업을 진행합니다.
            if os.path.exists(outcar_file):
#                os.system(f'echo "Merging OUTCAR in {dir}..."')
                merge_out(outcar_file)
                # 병합이 완료되면, gzip 명령어를 이용해 OUTCAR 파일을 다시 압축합니다.
                os.system(f"gzip {outcar_file}")
            
        """        
        for dir in dirs_err:
            outcar_file = os.path.join(dir, "OUTCAR")
            if os.path.exists(outcar_file):
                os.system(f'echo "Merging OUTCAR in {dir}..."')
                merge_out(outcar_file)
        """
        ########## OUTCAR Merge ######################

        """    
        # dirs = [d for d in os.listdir("./") if d.startswith("dir_POSCAR") and os.path.isdir(d)]
        dirs = sorted(
            [d for d in os.listdir("./") if d.startswith("dir_POSCAR") and os.path.isdir(d)],
            key=lambda x: int(x[10:])  # "con" 이후의 숫자로 정렬
        )
        """
        # "dir_POSCAR"로 시작하는 디렉토리 목록 가져오기
        dirs = [d for d in os.listdir("./") if d.startswith("dir_POSCAR") and os.path.isdir(d)]

        # 디렉토리가 2개 이상일 때만 정렬
        if len(dirs) > 1:
            dirs = sorted(dirs, key=lambda x: int(x[10:]))  # "dir_POSCAR" 뒤의 숫자로 정렬

        tot_dirs = len(dirs)
        first_check = True
        if idx == 1 and not os.path.exists("out.cfg"):
            open("out.cfg", "w").close()

        for idx_sub, d in enumerate(dirs, start=1):
            sub_path = os.path.join(d_step, d)
            print_progress_bar(idx_sub, tot_dirs, sub_path)
            pos_name = d.replace("dir_", "")

            outcar_file = os.path.join(d, "OUTCAR")
            outcar_gz   = outcar_file + ".gz"
            if os.path.exists(outcar_gz):
                os.system(f'gunzip {outcar_gz}')
            if not os.path.exists(outcar_file):
                con_index += 1
                first_check = False
                continue

            # POSCAR 6번째 줄 원소 정보 추출
            n_ele     = os.popen('sed -n 6p %s/POSCAR' % d).read().split()
            n_ele_pos = os.popen('sed -n 6p %s' % pos_name).read().split()

            # 첫 번째 실행 시 out.cfg 초기화
            if first_check and len(n_ele) < len(elist) and len(n_ele_pos) == len(elist) and idx == 1:
                open("out.cfg", "w").close()

            if len(n_ele) < len(elist) and len(n_ele_pos) == len(elist):
                os.system('echo " %s : %s / Modify out.cfg"' % (pos_name, " ".join(n_ele_pos)))
                shutil.copy("out.cfg", "out.cfg.prev")
                try:
                    os.system('%s/mlp convert-cfg %s/OUTCAR out.cfg --input-format=vasp-outcar --append >& /dev/null' % (mlip_path, d))
#                    os.system('%s/mlp convert-cfg %s/OUTCAR out.cfg --input-format=vasp-outcar --append --last >& /dev/null' % (mlip_path, d))  # Only last
                    update_feature_lines("out.cfg", root, con_values, dirs, con_index)
                    shutil.copy("out.cfg", "out.cfg.bak")
                    process_cfg_files("out.cfg.prev", "out.cfg", "out.cfg")
                    os.system(f'gzip {d}/OUTCAR')
                except subprocess.CalledProcessError as e:
                    print(f"\nError processing {d}: {e}")
            elif len(n_ele) < len(elist) and len(n_ele_pos) < len(elist):
                try:
                    os.system('%s/mlp convert-cfg %s/OUTCAR out.cfg --input-format=vasp-outcar --append >& /dev/null' % (mlip_path, d))
#                    os.system('%s/mlp convert-cfg %s/OUTCAR out.cfg --input-format=vasp-outcar --append --last >& /dev/null' % (mlip_path, d))  # Only last
                    update_feature_lines("out.cfg", root, con_values, dirs, con_index)
                    os.system(f'gzip {d}/OUTCAR')
                    os.system('echo " %s : %s / Keep out.cfg"' % (pos_name, " ".join(n_ele_pos)))
                except subprocess.CalledProcessError as e:
                    print(f"\nError processing {d}: {e}")
            else:
                try:
                    os.system('%s/mlp convert-cfg %s/OUTCAR out.cfg --input-format=vasp-outcar --append >& /dev/null' % (mlip_path, d))
#                    os.system('%s/mlp convert-cfg %s/OUTCAR out.cfg --input-format=vasp-outcar --append --last >& /dev/null' % (mlip_path, d))  # Only last
                    update_feature_lines("out.cfg", root, con_values, dirs, con_index)
                    os.system(f'gzip {d}/OUTCAR')
                except subprocess.CalledProcessError as e:
                    print(f"\nError processing {d}: {e}")

            con_index += 1
            first_check = False  

        ########## OUTCAR Restore ######################
        # 현재 디렉토리에서 error* 디렉토리 찾기 (error*.gz 파일이 있는 디렉토리만 선택)
        dirs_err = [
            dir for dir in sorted(glob.glob("dir_POSCAR*"), key=os.path.getmtime)
            if glob.glob(os.path.join(dir, "error*.gz"))  # 해당 디렉토리에 error*.gz 파일이 있는지 확인
        ]
    
        # OUTCAR 복원
        for dir in dirs_err:
            gz_file = os.path.join(dir, "OUTCAR.gz")
            outcar_file = os.path.join(dir, "OUTCAR")
            
            # 먼저, dir/OUTCAR.gz 파일이 존재하면 gunzip을 이용해 압축 해제합니다.
            if os.path.exists(gz_file):
                os.system(f"gunzip {gz_file}")
            
            # 압축 해제 후, OUTCAR 파일이 존재하면 병합 작업을 진행합니다.
            if os.path.exists(outcar_file):
#                os.system(f'echo "Restoring OUTCAR in {dir}..."')
                restore_outcar(outcar_file)  # merge_out 실행 후 OUTCAR 복원
                # 병합이 완료되면, gzip 명령어를 이용해 OUTCAR 파일을 다시 압축합니다.
                os.system(f"gzip {outcar_file}")
            """    
            outcar_file = os.path.join(dir, "OUTCAR")
            if os.path.exists(outcar_file):
                os.system(f'echo "Restoring OUTCAR in {dir}..."')
                restore_outcar(outcar_file)  # merge_out 실행 후 OUTCAR 복원
            """
        ########## OUTCAR MErge ######################

        subprocess.run(["cp", "-f", "out.cfg", base_dir], check=False)
        subprocess.run("rm -f out.cfg*", shell=True, check=False)
        os.chdir(base_dir)


Step_D()
print("\n")







