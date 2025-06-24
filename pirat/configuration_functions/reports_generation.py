import base64

def save_report(saving_path, name, current_time, raw_path_name, version, list_of_files, optional_saving_path=None):
    lista = ["Initial_distribution", "Initial_distribution_z_score",
             "5-to-5_overlap",
             "5-to-5_overlap_z_score", "Length_distribution_of_sequences_displaying_10nts_5-to-5_overlaps",
             "Length_distribution_of_sequences_displaying_10nts_5-to-5_overlaps_z_score",
             "heatmap", "seqlogo_pingpong"]
    svg_bytes = []
    svg_base64 = []
    for i in lista:
        with open(f"{saving_path}/{i}.svg", 'rb') as f:
            svg_bytes.append(f.read())
    for i in svg_bytes:
        svg_base64.append(base64.b64encode(i).decode())
    html_report = f"""
    <!DOCTYPE html>
    <html>
      <head>
        <title>piRAT {name} report</title>
        <style>
          .image-container {{
            display: flex;
      justify-content: center;
      margin: auto;
          }}
          .image-container img {{
            width: 100%;
            height: auto;
      justify-content: center;
          }}
        </style>
      </head>
      <body> 
        <h1>piRAT {name} report</h1>
        <p>piRAT {version}</p>
        <p>{current_time}</p>
        <p>Input path: {raw_path_name}</p>
        <p>Analyzed dataset: {list_of_files}</p>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[0]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[1]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[2]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[3]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[4]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[5]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[6]}"></div> 
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[7]}"></div> 
        </div>
      </body>
    </html>
    """
    if optional_saving_path is not None:
        with open(f"{optional_saving_path}/{name}_report.html", "w") as outputfile:
            outputfile.write(html_report)
    else:
        with open(f"{saving_path}/{name}_report.html", "w") as outputfile:
            outputfile.write(html_report)
    return None


def save_report_annotate(path, name, current_time, raw_path_name, version, list_of_files):
    lista = ["Initial_distribution", "Initial_distribution_z_score", "5-to-5_overlap",
             "5-to-5_overlap_z_score", "Length_distribution_of_the_sequences_displaying_10nts_5-to-5_overlaps",
             "Length_distribution_of_the_sequences_displaying_10nts_5-to-5_overlaps_z_score",
             "heatmap", "venn", "seqlogo_pingpong"]
    svg_bytes = []
    svg_base64 = []
    for i in lista:
        with open(f"{path}plots/{i}.svg", 'rb') as f:
            svg_bytes.append(f.read())
    for i in svg_bytes:
        svg_base64.append(base64.b64encode(i).decode())
    html_report = f"""
    <!DOCTYPE html>
    <html>
      <head>
        <title>piRAT annotate report of {name}</title>
        <style>
          .image-container {{
            display: flex;
      justify-content: center;
      margin: auto;
          }}
          .image-container img {{
            width: 100%;
            height: auto;
      justify-content: center;
          }}
        </style>
      </head>
      <body> 
        <h1>piRAT annotate report of {name}</h1>
        <p>piRAT {version}</p>
        <p>{current_time}</p>
        <p>Input path: {raw_path_name}</p>
        <p>Analyzed dataset: {list_of_files}</p>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[0]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[1]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[2]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[3]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[4]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[5]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[6]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[7]}"></div> 
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[8]}"></div> 
        </div>
      </body>
    </html>
    """

    with open(f"{path}reports/{name}_annotate_report.html", "w") as outputfile:
        outputfile.write(html_report)
    return None

def save_report_annotate_wo_primary(path, name, current_time, raw_path_name, version):
    lista = ["Initial_distribution", "Initial_distribution_z_score", "5-to-5_overlap",
             "5-to-5_overlap_z_score", "Length_distribution_of_the_sequences_displaying_10nts_5-to-5_overlaps",
             "Length_distribution_of_the_sequences_displaying_10nts_5-to-5_overlaps_z_score",
             "heatmap", "heatmap_pairs", "seqlogo_pingpong"]
    svg_bytes = []
    svg_base64 = []
    for i in lista:
        with open(f"{path}plots/{i}.svg", 'rb') as f:
            svg_bytes.append(f.read())
    for i in svg_bytes:
        svg_base64.append(base64.b64encode(i).decode())
    html_report = f"""
    <!DOCTYPE html>
    <html>
      <head>
        <title>piRAT annotate report of {name}</title>
        <style>
          .image-container {{
            display: flex;
      justify-content: center;
      margin: auto;
          }}
          .image-container img {{
            width: 100%;
            height: auto;
      justify-content: center;
          }}
        </style>
      </head>
      <body> 
        <h1>piRAT annotate report of {name}</h1>
        <p>piRAT {version}</p>
        <p>{current_time}</p>
        <p>Input path: {raw_path_name}</p>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[0]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[1]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[2]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[3]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[4]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[5]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[6]}"></div> 
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[7]}"></div> 
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[8]}"></div> 
        </div>
      </body>
    </html>
    """

    with open(f"{path}reports/{name}_annotate_report.html", "w") as outputfile:
        outputfile.write(html_report)
    return None

def save_final_report(path_sec, path_prim, path, name, range, automatic, current_time, raw_path_name, version, list_of_files, params):
    lista = ["Initial_distribution", "Initial_distribution_z_score", "5-to-5_overlap",
             "5-to-5_overlap_z_score", "Length_distribution_of_the_sequences_displaying_10nts_5-to-5_overlaps",
             "Length_distribution_of_the_sequences_displaying_10nts_5-to-5_overlaps_z_score",
             "heatmap", "heatmap_pairs", "venn", "seqlogo_pingpong"]
    lista_2 = ["Cluster_length_distribution.svg", "Cluster_length_distribution_quality.svg", "seqlogo_clusters.svg"]
    svg_bytes = []
    svg_base64 = []
    for i in lista:
        with open(f"{path_sec}plots/{i}.svg", 'rb') as f:
            svg_bytes.append(f.read())
    for i in lista_2:
        with open(f"{path_prim}plots/{i}", 'rb') as f:
            svg_bytes.append(f.read())
    for i in svg_bytes:
        svg_base64.append(base64.b64encode(i).decode())

    with open(f"{path_prim}other_data/Statistics_of_found_clusters.html", 'r') as f:
        statistics_html = f.read()
    range_autom = 'automatically' if automatic else 'manually'
    html_report = f"""
    <!DOCTYPE html>
    <html>
      <head>
        <title>piRAT final report of {name}</title>
         <style>
          .image-container {{
            display: flex;
            justify-content: center;
            margin: auto;
          }}
          .image-container img {{
            width: 100%;
            height: auto;
            justify-content: center;
          }}
        </style>
      </head>
      <body> 
        <h1>piRAT final report of {name}</h1>
        <p>piRAT {version}</p>
        <p>{current_time}</p>
        <p>Input path: {raw_path_name}</p>
        <p>Analyzed dataset: {list_of_files}</p>
        <h2>Analysis of mapped reads</h2>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[0]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[1]}"></div>
        </div>
        <div>
          <p style="text-align: center;">The piRNA length range was <b>{range_autom}</b> set as <b>[{range[0]}, {range[1]}] nts</b></p>
        </div>
        <h2>piRNA cluster annotation</h2>
        <p>Clustering configuration:</p>
        <ul>
            <li>MinReads (k): {params[0]}</li>
            <li>Eps: {params[1]}</li>
            <li>Range of size of piRNAs: {params[4]}</li>
            <li>variation threshold: {params[3]}</li>
            <li>Threads: {params[2]}</li>
        </ul>
        <div class="image-container">
          {statistics_html}
        </div>
        <h3>Cluster length distribution</h3>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[10]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[11]}"></div> 
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[12]}"></div> 
        </div>
        <h2>Ping-Pong signatures</h2>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[2]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[3]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[4]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[5]}"></div>
        </div>
        <div>
          <p style="text-align: center;">Sequences with 10 nts 5'-to-5' overlap and within the <b>{range_autom}</b> defined piRNA range <b>[{range[0]}, {range[1]}] nts</b> are considered piRNAs</p>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[6]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[7]}"></div> 
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[8]}"></div> 
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[9]}"></div> 
        </div>
      </body>
    </html>
    """
    with open(f"{path}reports/{name}_final_report.html", "w") as outputfile:
        outputfile.write(html_report)
    with open(f"{path}{name}_final_report.html", "w") as outputfile:
        outputfile.write(html_report)
    return None


def save_final_report_clustering(path_prim, path, name, current_time, raw_path_name, version, list_of_files,  params):
    lista_2 = ["Cluster_length_distribution.svg", "Cluster_length_distribution_quality.svg", "seqlogo_clusters.svg"]
    svg_bytes = []
    svg_base64 = []
    for i in lista_2:
        with open(f"{path_prim}plots/{i}", 'rb') as f:
            svg_bytes.append(f.read())
    for i in svg_bytes:
        svg_base64.append(base64.b64encode(i).decode())

    with open(f"{path_prim}other_data/Statistics_of_found_clusters.html", 'r') as f:
        statistics_html = f.read()
    html_report = f"""
    <!DOCTYPE html>
    <html>
      <head>
        <title>piRAT final report of clustering of {name}</title>
         <style>
          .image-container {{
            display: flex;
            justify-content: center;
            margin: auto;
          }}
          .image-container img {{
            width: 100%;
            height: auto;
            justify-content: center;
          }}
        </style>
      </head>
      <body> 
        <h1>piRAT final report of clustering of {name}</h1>
        <p>piRAT {version}</p>
        <p>{current_time}</p>
        <p>Input path: {raw_path_name}</p>
        <p>Analyzed dataset: {list_of_files}</p>
        <p>Clustering configuration:</p>
        <ul>
            <li>MinReads (k): {params[0]}</li>
            <li>Eps: {params[1]}</li>
            <li>Range of size of piRNAs: {params[4]}</li>
            <li>variation threshold: {params[3]}</li>
            <li>Threads: {params[2]}</li>
        </ul>
        <div class="image-container">
          {statistics_html}
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[0]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[1]}"></div> 
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[2]}"></div> 
        </div>
      </body>
    </html>
    """

    with open(f"{path}reports/{name}_final_report_clustering.html", "w") as outputfile:
        outputfile.write(html_report)
    return None

def save_report_no_ping_pong(saving_path, name, current_time, raw_path_name, version):
    lista = ["Initial_distribution", "Initial_distribution_z_score",
             "5-to-5_overlap",
             "5-to-5_overlap_z_score"]
    svg_bytes = []
    svg_base64 = []
    for i in lista:
        with open(f"{saving_path}/{i}.svg", 'rb') as f:
            svg_bytes.append(f.read())
    for i in svg_bytes:
        svg_base64.append(base64.b64encode(i).decode())
    html_report = f"""
    <!DOCTYPE html>
    <html>
      <head>
        <title>{name} report</title>
        <style>
          .image-container {{
            display: flex;
      justify-content: center;
      margin: auto;
          }}
          .image-container img {{
            width: 100%;
            height: auto;
      justify-content: center;
          }}
        </style>
      </head>
      <body> 
        <h1>{name} report</h1>
        <p>piRAT {version}</p>
        <p>{current_time}</p>
        <p>Input path: {raw_path_name}</p>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[0]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[1]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[2]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[3]}"></div>
        </div>
      </body>
    </html>
    """

    with open(f"{saving_path}{name}_report.html", "w") as outputfile:
        outputfile.write(html_report)
    return None

def save_final_report_no_ping_pong(path_sec, path_prim, path, name, range, automatic, current_time, raw_path_name, version):
    lista = ["Initial_distribution", "Initial_distribution_z_score", "5-to-5_overlap",
             "5-to-5_overlap_z_score", "heatmap", "venn"]
    lista_2 = ["Cluster_length_distribution.svg", "Cluster_length_distribution_quality.svg", "seqlogo_clusters.svg"]
    svg_bytes = []
    svg_base64 = []
    for i in lista:
        with open(f"{path_sec}plots/{i}.svg", 'rb') as f:
            svg_bytes.append(f.read())
    for i in lista_2:
        with open(f"{path_prim}plots/{i}", 'rb') as f:
            svg_bytes.append(f.read())
    for i in svg_bytes:
        svg_base64.append(base64.b64encode(i).decode())

    with open(f"{path_prim}other_data/Statistics_of_found_clusters.html", 'r') as f:
        statistics_html = f.read()
    range_autom = 'automatically' if automatic else 'manually'
    html_report = f"""
    <!DOCTYPE html>
    <html>
      <head>
        <title>piRAT final  report of {name}</title>
         <style>
          .image-container {{
            display: flex;
            justify-content: center;
            margin: auto;
          }}
          .image-container img {{
            width: 100%;
            height: auto;
            justify-content: center;
          }}
        </style>
      </head>
      <body> 
        <h1>piRAT final  report of {name}</h1>
        <p>piRAT {version}</p>
        <p>{current_time}</p>
        <p>Input path: {raw_path_name}</p>
        <h2>Analysis of mapped reads</h2>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[0]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[1]}"></div>
        </div>
        <div>
          <p style="text-align: center;">The piRNA length range was <b>{range_autom}</b> set as <b>[{range[0]}, {range[1]}] nts</b></p>
        </div>
        <h2>piRNA cluster annotation</h2>
        <div class="image-container">
          {statistics_html}
        </div>
        <h3>Cluster length distribution</h3>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[6]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[7]}"></div> 
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[8]}"></div> 
        </div>
        <h2>Ping-Pong signatures</h2>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[2]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[3]}"></div>
        </div>
        <div>
          <p style="text-align: center;">Sequences with 10 nts 5'-to-5' overlap and within the <b>{range_autom}</b> defined piRNA range <b>[{range[0]}, {range[1]}] nts</b> are considered piRNAs</p>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[4]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[5]}"></div> 
        </div>
      </body>
    </html>
    """

    with open(f"{path}reports/{name}_final_report.html", "w") as outputfile:
        outputfile.write(html_report)
    return None

def save_report_annotate_wo_primary_and_secondary(path, name, current_time, raw_path_name, version):
    lista = ["Initial_distribution", "Initial_distribution_z_score", "5-to-5_overlap",
             "5-to-5_overlap_z_score", "heatmap", "heatmap_pairs"]
    svg_bytes = []
    svg_base64 = []
    for i in lista:
        with open(f"{path}plots/{i}.svg", 'rb') as f:
            svg_bytes.append(f.read())
    for i in svg_bytes:
        svg_base64.append(base64.b64encode(i).decode())
    html_report = f"""
    <!DOCTYPE html>
    <html>
      <head>
        <title>piRAT annotate report of {name}</title>
        <style>
          .image-container {{
            display: flex;
      justify-content: center;
      margin: auto;
          }}
          .image-container img {{
            width: 100%;
            height: auto;
      justify-content: center;
          }}
        </style>
      </head>
      <body> 
        <h1>piRAT annotate report of {name}</h1>
        <p>piRAT {version}</p>
        <p>{current_time}</p>
        <p>Input path: {raw_path_name}</p>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[0]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[1]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[2]}"></div> 
           <div><img src="data:image/svg+xml;base64,{svg_base64[3]}"></div>
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[4]}"></div> 
        </div>
        <div class="image-container">
           <div><img src="data:image/svg+xml;base64,{svg_base64[5]}"></div> 
        </div>
      </body>
    </html>
    """

    with open(f"{path}reports/{name}_annotate_report.html", "w") as outputfile:
        outputfile.write(html_report)
    return None