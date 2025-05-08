import subprocess
import os
import pandas as pd

# Definir o diretório onde o MadGraph está instalado e o diretório de trabalho
madgraph_path = "/mnt/d/codigos/posdoc/resonant-scalar-madg/madgraph/MG5_aMC_v3_6_0/bin/mg5_aMC"
working_dir = "/mnt/d/codigos/posdoc/resonant-scalar-madg/madgraph/workdir"

# Crie um arquivo de comando do MadGraph para geração de eventos
command_file = os.path.join(working_dir, "process_command.txt")
with open(command_file, "w") as f:
    f.write("""
    generate e+ e- > mu+ mu-
    output generated_events
    launch
    set run_card detector on    # Ativa a conversão para ROOT
    """)
    
# Execute o MadGraph usando o subprocess
try:
    subprocess.run([madgraph_path, command_file], cwd=working_dir, check=True)
    print("MadGraph executado com sucesso.")
except subprocess.CalledProcessError as e:
    print("Erro ao executar MadGraph:", e)

# Caminho para o arquivo de saída com eventos gerados (lhef ou root, dependendo da configuração)
output_file = os.path.join(working_dir, "generated_events/Events/run_01/unweighted_events.lhe")
output_zip_file = os.path.join(working_dir, "generated_events/Events/run_01/unweighted_events.lhe.gz")

# Função para ler o arquivo LHE e converter para um DataFrame do Pandas
def parse_lhe_file(file_path):
    events = []
    subprocess.run(["gzip", "-d", output_zip_file], cwd=working_dir)
    print("Arquivo LHE descompactado com sucesso.")
    with open(file_path, 'r') as file:
        in_event = False
        for line in file:
            if "<event>" in line:
                in_event = True
                event = []
            elif "</event>" in line:
                in_event = False
                events.append(event)
            elif in_event:
                data = line.strip().split()
                if len(data) >= 6:
                    particle_data = list(map(float, data))
                    event.append(particle_data)
    # Converte a lista de eventos para um DataFrame
    columns = ["PID", "Status", "Mother1", "Mother2", "Color1", "Color2", "Px", "Py", "Pz", "E", "M", "Lifetime", "Spin"]
    df = pd.DataFrame([item for sublist in events for item in sublist], columns=columns)
    return df

# Parse o arquivo LHE
df_events = parse_lhe_file(output_file)

# Mostrar os primeiros eventos
print(df_events.head())
