import os
import ssl
from time import sleep, time
from Bio import Entrez
import requests
from requests.exceptions import HTTPError
import numpy as np
import matplotlib.pyplot as plt
import json

# Security key to use the database NCBI

if (not os.environ.get('PYTHONHTTPSVERIFY', '') and
        getattr(ssl, '_create_unverified_context', None)):
    ssl._create_default_https_context = ssl._create_unverified_context


# ---------------------------------------------------------------------------------
#                   Define the Input and Output Files
# ---------------------------------------------------------------------------------

def define_files():
    print('=-' * 15)
    print('{:^30}'.format('Welcome on Board Sailor !'))
    print('-=' * 15)
    # Parameters to access the NCBI service:
    Entrez.tool = 'SHIP'  # Tool that is accessing
    Entrez.email = ' '  # Tool user email - Always tell NCBI who you are.
    while True:
        Entrez.email = str(input('Enter your email to access NCBI: '))
        if '@' not in Entrez.email:
            print('\033[1:31mInvalid E-mail\033[m')
        else:
            break

    # GFF file entry for analysis
    while True:
        name_file_read = str(input('Enter the path and file name: ')).strip()
        if name_file_read[-3:] != 'gff' and name_file_read[-4:] != 'gff3':  # Checks whether the input file is in gff
            print('Inform file in gff format')
        else:
            try:
                file_r = open(name_file_read, 'r')  # Function opens a file, and returns it as a file object.
            except FileNotFoundError:
                print('File not found, please check the path or file name with the extension.')
            except Exception as error:
                print(error)
            else:
                break

    interval = input('Enter the interval for genomic analysis (default 500): ')

    if len(interval) == 0:
        return file_r, name_file_read, 500
    else:
        return file_r, name_file_read, int(interval)


# ---------------------------------------------------------------------------------
#                    PARSE of the GFF File
# ---------------------------------------------------------------------------------

def parse_gff(file_r):
    print('--' * 15)
    print('ANALYZING ...')
    print('--' * 15)
    print('''Performing the following steps:
    [1] Selecting and sorting the genes by their start within each chromosome.
    [2] Disregarding completely overlapping genes.
    [3] Filtering out convergent, divergent and tandem genes.''')
    print('--' * 15)

    file_format = 1  # 1 indicates the standard model and 2 indicates Ensembl model
    cont_line = 0
    result = list()
    cont_chromosome = 0
    analyzed_type = {}  # Stores the types to be analyzed
    unanalyzed_type = {}  # Stores the types that will not be analyzed

    with open('features.json', 'r') as f:
        data_json = json.load(f)
    print(f'features : {list(data_json.keys())}')
    print('--' * 15)

    #for line in file_r:
    lista_arquivo = file_r.readlines() # Transforma o arquivo em uma Lsita
    chrom_anterior = ''
    for pos_i in range (0,len(lista_arquivo)):
        line = lista_arquivo[pos_i]
        li = line.split('\t')  # Transforms a TAB-separated line into a LIST for easy manipulation

        if pos_i + 1 < len(lista_arquivo): # pega a proxima linha da lista
            lp = lista_arquivo[pos_i+1].split('\t')

        if 'Alias=' in line:  # Identifying the chromosome and renaming column 0
            file_format = 2  # Ensembl model
            alias = li[8]
            start_position = alias.index("Alias=") + 6
            end_position = alias.find(',', start_position)  # If the chromosome id is between Alias= and ','
            if end_position == -1:
                end_position = alias.find(';', start_position)  # If the chromosome id is between Alias= and ';'
            if end_position == -1:
                end_position = alias.find('\t')  # If the chromosome id is between Alias= and end of the line

            id_chromosome = alias[start_position:end_position]  # Get the chromosome id
            print(f'{id_chromosome}, ', end='')  # Shows the chromosome that is being analyzed
            cont_chromosome += 1
            if cont_chromosome == 5:
                print('')
                cont_chromosome = 0
        find_feature = False
        if line[0] != '#':
            feature = li[2]
            for k, v in data_json.items():
                if k in li[2] and len(k) == len(li[2]):
                    if type(v) is list and lp[2] in v:
                        find_feature = True
                    else:

                        find_feature = True
                    break

            if find_feature:
                if feature in analyzed_type:
                    analyzed_type[feature] += 1
                else:
                    analyzed_type[feature] = 1

                if file_format == 2:  # Ensemble model
                    li[0] = id_chromosome
                    li.append(line[0])
                else:
                    li.append('N')  # File with a different format than Ensemble

                id_chromosome = li[0]
                if id_chromosome != chrom_anterior:
                    print(f'{li[0]} ', end='')
                    chrom_anterior = id_chromosome


                li[3] = int(li[3])  # Save the start and end columns of the genes as integers
                li[4] = int(li[4])  # To perform the calculations
                #  Insert line
                result.append(li[:])
                li.clear()
            else:
                if feature in unanalyzed_type:
                    unanalyzed_type[feature] += 1
                else:
                    unanalyzed_type[feature] = 1
            cont_line += 1

    print('\n')
    print(f'{" Types that will be processed ":-^50}')
    for k, q in analyzed_type.items():
        print(f'{k} = {q}')
    print(f'{" Types that will not be processed ":-^50}')
    for k, q in unanalyzed_type.items():
        print(f'{k} = {q}')
    print('-' * 50)
    file_r.close()
    return result


# ---------------------------------------------------------------------------------
#                    Remove Total Overlap
# ---------------------------------------------------------------------------------

def remove_overlap(result):
    def fucKey(e):
        return e[0], e[3]

    result.sort(key=fucKey)  # Sort genes by chromosome number and Start
    result_without_overlap = []  # Looking for and removing total overlap
    cl = 0  # Current line
    current_chromosome = result[cl][0]

    try:
        for l1 in range(1, len(result)):
            posterior_chromosome = result[l1][0]
            if current_chromosome == posterior_chromosome:  # Check if the genes are on the same chromosome
                if result[l1][3] >= result[cl][3] and result[l1][4] <= result[cl][4]:  # Finding overlays
                    pass
                else:
                    cl = l1
                    result_without_overlap.append(result[l1])  # Add genes that are not totally overlapping
            else:
                current_chromosome = posterior_chromosome
                cl = l1

    except Exception as error:
        print(f'Overlap error = {error}')

    amount_overlap = len(result) - len(result_without_overlap)
    return result_without_overlap, amount_overlap


# ---------------------------------------------------------------------------------
#                    Assemble Table
# ---------------------------------------------------------------------------------

def assemble_table(non_overlap_list, interval, amount_genes, amount_overlap):
    import collections as coll
    amount_chromosome = 1
    end_list = len(non_overlap_list) - 1
    ranges_dictionary = {
        0: {'description': 'Overlapping',
            'amount_tandem': 0,
            'amount_divergent': 0,
            'amount_convergent': 0}}
    for l1 in range(0, end_list):
        l2 = l1 + 1
        if l2 < end_list:
            if non_overlap_list[l1][0] != non_overlap_list[l2][0]:  # Check for chromosome change
                l1 = l2
                l2 = l1 + 1
                amount_chromosome += 1
            # Calculates the difference between the start of the second gene and the end of the first gene

            difference = non_overlap_list[l2][3] - non_overlap_list[l1][4]
            if difference > 0:  # Checks for gaps between genes
                quotient = difference // interval
                rest = difference % interval
                if rest > 0:
                    key = quotient + 1
                    if key in ranges_dictionary:  # Checks if the range key is already in the dictionary.
                        position = ranges_dictionary.get(key)
                        if (non_overlap_list[l1][6] == '-' and non_overlap_list[l2][6] == '-') or (
                                non_overlap_list[l1][6] == '+' and non_overlap_list[l2][6] == '+'):  # if it's tandem
                            quantity = position.get('amount_tandem')
                            quantity += 1
                            position.update({'amount_tandem': quantity})
                        if non_overlap_list[l1][6] == '-' and non_overlap_list[l2][6] == '+':  # If it's divergent
                            quantity = position.get('amount_divergent')
                            quantity += 1
                            position.update({'amount_divergent': quantity})
                        if non_overlap_list[l1][6] == '+' and non_overlap_list[l2][6] == '-':  # If it's convergent
                            quantity = position.get('amount_convergent')
                            quantity += 1
                            position.update({'amount_convergent': quantity})
                    else:  # Adding a new element (interval) to the dictionary
                        start_value_description = interval * quotient + 1
                        end_value_description = start_value_description + interval - 1
                        description = f'{float(start_value_description):,.0f}-{float(end_value_description):,.0f} bp'
                        if (non_overlap_list[l1][6] == '-' and non_overlap_list[l2][6] == '-') or (
                                non_overlap_list[l1][6] == '+' and non_overlap_list[l2][6] == '+'):  # if it's tandem
                            ranges_dictionary[key] = {'description': description, 'amount_tandem': 1,
                                                      'amount_divergent': 0, 'amount_convergent': 0}
                        if non_overlap_list[l1][6] == '-' and non_overlap_list[l2][6] == '+':  # If it's divergent
                            ranges_dictionary[key] = {'description': description, 'amount_tandem': 0,
                                                      'amount_divergent': 1, 'amount_convergent': 0}
                        if non_overlap_list[l1][6] == '+' and non_overlap_list[l2][6] == '-':  # If it's convergent
                            ranges_dictionary[key] = {'description': description, 'amount_tandem': 0,
                                                      'amount_divergent': 0, 'amount_convergent': 1}
            else:  # difference <= 0: Overlapping
                position = ranges_dictionary.get(0)
                if (non_overlap_list[l1][6] == '-' and non_overlap_list[l2][6] == '-') or (
                        non_overlap_list[l1][6] == '+' and non_overlap_list[l2][6] == '+'):  # if it's tandem
                    quantity = position.get('amount_tandem')
                    quantity += 1
                    position.update({'amount_tandem': quantity})
                if non_overlap_list[l1][6] == '-' and non_overlap_list[l2][6] == '+':  # If it's divergent
                    quantity = position.get('amount_divergent')
                    quantity += 1
                    position.update({'amount_divergent': quantity})
                if non_overlap_list[l1][6] == '+' and non_overlap_list[l2][6] == '-':  # If it's convergent
                    quantity = position.get('amount_convergent')
                    quantity += 1
                    position.update({'amount_convergent': quantity})

    ordered_dictionary = coll.OrderedDict(sorted(ranges_dictionary.items()))
    contTandem = contDiverge = contConverge = 0

    print(f'{"Interval (bp)":^20} {"Tandem":>20} {"Divergent":>20} {"Convergent":>20}')
    for k, v in ordered_dictionary.items():
        for key, value in v.items():
            print(f'{value:>20}', end='')
            if key == 'amount_tandem':
                contTandem += value
            if key == 'amount_divergent':
                contDiverge += value
            if key == 'amount_convergent':
                contConverge += value
        print()

    '''
    print()
    print(f'Quantity of Chromosomes = {amount_chromosome}')
    print(f'Number of genes located = {amount_genes}')
    print(f'Total intergenic intervals = {contDiverge + contConverge + contTandem}')
    print(f'Total Intervals flanked by Tandem genes     = {contTandem}')
    print(f'Total Intervals flanked by Divergent genes  = {contDiverge}')
    print(f'Total Intervals flanked by Convergent genes = {contConverge}')
    print(f'Number of complete overlaps between genes   = {amount_overlap}\n')
    print()
    '''
    print(f'''
    Quantity of Chromosomes = {amount_chromosome} 
    Number of genes located = {amount_genes}
    Total intergenic intervals = {contDiverge + contConverge + contTandem}
    Total Intervals flanked by Tandem genes     = {contTandem}
    Total Intervals flanked by Divergent genes  = {contDiverge}
    Total Intervals flanked by Convergent genes = {contConverge}
    Number of complete overlaps between genes   = {amount_overlap}
    
    ''')
    return ordered_dictionary


# ----------------------------------------------------
#     Plotagem dos dados
# ---------------------------------------------------
def plot_genes(dicto_genes, interval_plot, file_name):
    soma_tamdem = soma_divergente = soma_convergente = 0
    dict_final = dict()
    labels = list()
    list_results = list()
    cont = 0  # conta o nivel a ser plotado
    # monta um novo dicionario incluindo somente os itens menores que o interval_plot
    # e inclui um novo item com a somatoria dos itens maiores que o interval_plot
    for k, v in dicto_genes.items():
        if interval_plot <= cont:
            for ki, vi in v.items():
                if k == interval_plot:
                    if ki == 'description':
                        fim = vi.index("-")
                        descricao = '≥'
                        descricao += vi[0:fim]
                        descricao += ' bp'
                if ki == 'amount_tandem':
                    soma_tamdem += vi
                elif ki == 'amount_divergent':
                    soma_divergente += vi
                elif ki == 'amount_convergent':
                    soma_convergente += vi
        else:
            dict_final[k] = v
        cont += 1

    dict_final[interval_plot] = {'description': descricao, 'amount_tandem': soma_tamdem,
                                 'amount_divergent': soma_divergente, 'amount_convergent': soma_convergente}

    # aqui começa a plotagem
    list_dados = dict_final.values()
    results = list()
    # dados = np.array(list(dados.values()))
    for v in list_dados:
        for ki, vi in v.items():
            if ki == 'description':
                labels.append(vi)
            else:
                list_results.append(vi)
        results.append(list_results[:])
        list_results.clear()

    # labels = list(results.keys())
    data = np.array(results)
    category_names = ['Tandem', 'Divergent', 'Convergent']
    data_cum = data.cumsum(axis=1)
    # category_colors = plt.get_cmap('RdYlBu')(np.linspace(0.15, 0.85, data.shape[1]))

    category_colors = plt.get_cmap('RdYlBu')(np.linspace(0.20, 0.80, data.shape[1]))
    # category_colors = ['r', 'g', 'b']
    # Cria  a figura e o Axes para a plotagem
    # fig, ax = plt.subplots(figsize=(15,10))
    fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(18, 6))
    lim_x = np.sum(data, axis=1).max() / 2 + 1000
    fig.suptitle('', fontsize=16)
    SMALL_SIZE = 11
    MEDIUM_SIZE = 15
    BIGGER_SIZE = 25

    plt.rc('font', size=MEDIUM_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels (Numbers y)
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels (Numbers y)
    plt.rc('legend', fontsize=BIGGER_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)
    colors_axes = ['royalblue', 'forestgreen', 'crimson']
    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        ax[i].invert_yaxis()
        ax[i].xaxis.set_visible(True)
        ax[i].set_xlim(0, lim_x)
        if i > 0:
            ax[i].yaxis.set_visible(False)
        widths = data[:, i]
        starts = data_cum[:, i] - widths

        ax[i].barh(labels, widths, left=0, height=0.5, label=colname, color=colors_axes[i])
        labels_x = ax[i].get_xticklabels()
        plt.setp(labels_x, rotation=45, horizontalalignment='right')
        xcenters = widths / 2
        r, g, b, _ = color
        # text_color = 'white' if r * g * b < 0.8 else 'darkgrey'
        for y, (x, c) in enumerate(zip(xcenters, widths)):
            ax[i].text(x, y, f'{c:,.0f}', ha='left', va='center', color='k')

        ax[i].legend(ncol=len(category_names), bbox_to_anchor=(0, 1),

                     loc='lower left', fontsize='small')
    file_name = file_name[:-4]
    file_name += '_plot.png'
    plt.savefig(file_name)
    # plt.show()
    # plt.close('all')  # all open plots are correctly closed after each run


# ---------------------------------------------------------------------------------
#                    Generate the FASTA of the Selected Intervals
# ---------------------------------------------------------------------------------

def fasta(non_overlap_list, file_name):
    regulatory_list = []
    regulatory_dict = dict()
    list_chromosome = list()
    choice_Cross_References : str = ''
    choice_regulatory: str = ''
    species: str = ''

    name_file_write = ''  # Preparing to create the writing file:
    try:
        name_file_write = file_name[:-4]
        name_file_write += '_result_SHIP.txt'  # Mount the output file name by changing the ".gff" to "_result_SHIP.txt"
        file_w = open(name_file_write, 'w')
    except Exception as error:
        print(error)

    print('--' * 15)
    print('Obtaining the FASTA Sequence')
    print('--' * 15)

    # Specification of the type of genes that flank the intergenic region:
    first_gene = ''
    second_gene = ''

    while True:
        print('Types of flanking genes: \n[1] CONVERGENT (+ -) \n[2] DIVERGENT (- +) \n[3] TANDEM (+ +/- -)')
        choice = int(input('> Enter your choice: '))
        try:
            if choice not in (1, 2, 3):
                print('Choose 1, 2 or 3')
            if choice == 1:
                first_gene = '+'
                second_gene = '-'
                break
            if choice == 2:
                first_gene = '-'
                second_gene = '+'
                break
            if choice == 3:
                break
        except Exception as error:
            print(error)
            print('Just type numbers!')

    # Specification of the range size of the intergenic region:
    while True:
        value1 = str(input('> Minimum size of the intergenic region (Standard 1500): '))
        value2 = str(input('> Maximum size of the intergenic region (Standard 2000): '))
        try:
            if value1 == '':
                value1 = 0
            if value2 == '':
                value2 = 0
            ranger_min = int(value1)
            ranger_max = int(value2)
            if ranger_max == 0 and ranger_min == 0:
                ranger_max = 2000
                ranger_min = 1500
                break
            elif ranger_max <= ranger_min:
                print('The minimum value must be less than the maximum value.')
            else:
                break
        except Exception as error:
            print(error)
            print('Just type numbers!')

    # -------------------------------------------------------------------
    #    Define whether California data analysis will be performed
    # -------------------------------------------------------------------

    while True:
        california_analysis = input('Do you want to perform Data Base California analysis ? (Y/N) ').strip().upper()[0]
        if california_analysis not in ("Y", "N"):
            print('Choose Y or N')
        else:
            break

    # -------------------------------------------------------------------
    #    Define whether regulation data analysis will be performed
    # -------------------------------------------------------------------

    while True:
        choice_Cross_References = input('Do you want to perform Cross References? (Y/N) ').strip().upper()[0]
        ## verficar se foi digitado Y ou N

        choice_regulatory = input('Do you want to perform regulatory analysis? (Y/N) ').strip().upper()[0]
        if choice_regulatory == 'Y':
            species = str(input('Enter the species to be analyzed: '))
            read_file_name = str(input('Enter the regulatory build Path and File name: ')).strip()
            # Checks whether the input file is in gff format
            if read_file_name[-3:] != 'gff' and read_file_name[-4:] != 'gff3':
                print('Inform file in gff format')
            else:
                try:
                    file_r = open(read_file_name)
                except FileNotFoundError:
                    print('File not found, please check the path or file name with the extension.')
                except Exception as erro:
                    print(erro)
                else:
                    for line in file_r:  # Copying the file into a list
                        regulatory_list.append(line.split('\t'))

                    def fucKey(e):
                        return e[0], e[3]

                    regulatory_list.sort(key=fucKey)  # Sort genes by chromosome number and Start

                    current_chromosome = regulatory_list[0][0]
                    # Creating a dictionary where the key will be the chromosome id and the element a list of all
                    # chromosomes related to the key.
                    for line in regulatory_list:
                        if line[0] == current_chromosome:
                            list_chromosome.append(line[:])
                        else:
                            regulatory_dict[current_chromosome] = list_chromosome[:]
                            list_chromosome.clear()
                            current_chromosome = line[0]
                    break
        else:
            break

    print('--' * 15)
    print('PROCESSING ...')
    print('--' * 15)
    print(f'[1] Selecting the intergenic regions with size between {ranger_min} to {ranger_max} bp.')
    print('[2] Searching for sequences at NCBI.')
    print('[3] Crossing information with other databases.')
    print('Wait for processing... Go have a cup of coffee.')
    print('⤍', end=' ')
    cell_types_dict = dict()  # Used to build a dictionary by grouping cell types by activity
    cont_interval = 0

    ## tem colocar a verificação para testar se o arquivo realmente foi lido
    with open('tracks.json', 'r') as f:
        track_json = json.load(f)


    try:
        for l1 in range(0, len(non_overlap_list)):
            process = False
            l2 = l1 + 1
            if l2 < len(non_overlap_list):  # Check if the types of flanking genes are the ones the user has chosen
                if choice == 3:
                    if (non_overlap_list[l1][6] == '+' and non_overlap_list[l2][6] == '+') or (
                            non_overlap_list[l1][6] == '-' and non_overlap_list[l2][6] == '-'):
                        process = True
                else:
                    if non_overlap_list[l1][6] == first_gene and non_overlap_list[l2][6] == second_gene:
                        process = True

                if process:
                    difference = non_overlap_list[l2][3] - non_overlap_list[l1][4]
                    if ranger_min <= difference <= ranger_max:  # Selecting the size of the chosen range
                        try:
                            # v == value
                            v_id = non_overlap_list[l1][0]
                            v_seq_start = str(int(non_overlap_list[l1][4]) + 1)
                            v_seq_stop = str(int(non_overlap_list[l2][3]) - 1)
                            # Accessing the NCBI database through the Entrez API
                            handle = Entrez.efetch(db='nuccore', id=v_id, seq_start=v_seq_start, seq_stop=v_seq_stop,
                                                   rettype='fasta', retmode='xml')
                            records = Entrez.parse(handle)

                            # Assembling the FASTA format.
                            for record in records:
                                line_cromossomo = record.get('TSeq_defline')
                                ## pegando o chromosso par pesquisa na base de dados da california
                                pos_chrom_ini = line_cromossomo.index("chromosome") + len("chromosome") + 1
                                pos_chrom_end = line_cromossomo.find(',', pos_chrom_ini)
                                num_chrom = 'chr' + str(line_cromossomo[pos_chrom_ini:pos_chrom_end]).strip()
                                # print(f'num_chrom = {num_chrom}')
                                ## aqui começa a verificação na Base de dados da California
                                print_fast = False
                                if california_analysis == 'Y':
                                    for k_track, v_track in track_json.items():
                                        print(f'k_track = {k_track}')
                                        if k_track == 'genome':
                                            genome_search = v_track
                                        if k_track != 'genome':
                                            track = k_track
                                            server = "https://api.genome.ucsc.edu/"
                                            ext = f"/getData/track?genome={genome_search};chrom={num_chrom};track={track};start={v_seq_start};end={v_seq_stop};maxItemsOutput=-1"
                                            # https://api.genome.ucsc.edu/getData/track?genome=sacCer3;track=sgdOther;chrom=chrXII;start=150164;end=151388;maxItemsOutput=-1
                                            print(server + ext)
                                            r = requests.get(server + ext, headers={"Content-Type": "application/json"})
                                            r.raise_for_status()
                                            jsonResponse = r.json()
                                            print(jsonResponse)
                                            pega_trilha = jsonResponse.get(track)
                                            # print(f'pega_trilha = {pega_trilha}')
                                            print(f'v_track = {v_track}')
                                            if pega_trilha != None:
                                                # se o retorno o BD da california vier vazio
                                                if len(pega_trilha) == 0:
                                                    print_fast = True
                                                # se o retorno o BD da california nao vier vazio mais a trilha no json for vazia
                                                elif len(v_track) == 0 and len(pega_trilha) > 0:
                                                    print_fast = False

                                                elif len(v_track) > 0 and len(pega_trilha) > 0:
                                                    for val in v_track:
                                                        if str(val).upper().strip() in str(pega_trilha).upper().strip():
                                                            # print(f'teve casamento aqui ... :{str(val).upper().strip()}')
                                                            print_fast = False
                                                            break
                                            if print_fast == False:
                                                break


                                # print(f'california_analysis = {california_analysis} Print_fasta = {print_fast} ')
                                if print_fast:
                                    line = ''
                                    file_w.write(
                                        '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
                                    file_w.write('\n')
                                    file_w.write(str(non_overlap_list[l1]))
                                    file_w.write('\n')
                                    file_w.write(str(non_overlap_list[l2]))
                                    file_w.write('\n')
                                    file_w.write(
                                        '------------------------------------------------------------------------------------')
                                    file_w.write('\n')
                                    file_w.write('Length = ' + str(difference - 1) + 'bp')
                                    file_w.write('\n')
                                    line = '>' + record.get('TSeq_accver')
                                    line += ':'
                                    line += v_seq_start
                                    line += '-'
                                    line += v_seq_stop
                                    line += ' '
                                    line += record.get('TSeq_defline')
                                    file_w.write(line)
                                    file_w.write('\n')
                                    line_fasta = record.get('TSeq_sequence')
                                    file_w.write(line_fasta)
                                    file_w.write('\n')
                                    print('☕︎︎', end='')

                                    sleep(5)  # 2 second delay to avoid overloading the server
                                    if choice_regulatory == 'Y':
                                        # ------------------------------------------------------------------
                                        #  Do research on Ensembl
                                        # -------------------------------------------------------------------
                                        start_range = non_overlap_list[l1][3]
                                        end_range = non_overlap_list[l2][4]
                                        if non_overlap_list[9] != 'N':
                                            id_chromosome = non_overlap_list[l1][9]  # Ensemble format file
                                        else:
                                            id_chromosome = non_overlap_list[l1][0]  # File in standard format

                                        list_key_chromosome = regulatory_dict[id_chromosome]
                                        file_w.write('\n')
                                        file_w.write(
                                            '------------------------------------------------------------------------------------')
                                        file_w.write('\n')
                                        file_w.write('Regulation: ')
                                        file_w.write('\n')
                                        file_w.write(
                                            '------------------------------------------------------------------------------------')
                                        file_w.write('\n')

                                        for line in list_key_chromosome:
                                            if start_range <= int(line[3]) <= end_range:
                                                # Looking for id
                                                attributes = line[8]
                                                start_position = attributes.index(":") + 1
                                                end_position = attributes.find(';', start_position)
                                                id_search = attributes[start_position:end_position]
                                                # Mounting query URL in the Ensemble Regulation
                                                server = "https://rest.ensembl.org"
                                                ext = f"/regulatory/species/{species}/id/{id_search}?activity=1"
                                                try:
                                                    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
                                                    print(server + ext)
                                                    r.raise_for_status()
                                                    jsonResponse = r.json()
                                                    # file_w.write(str(jsonResponse))
                                                    file_w.write('\n')
                                                    file_w.write(str(line))
                                                    file_w.write('\n')
                                                    for value in jsonResponse:
                                                        if type(value) is dict:
                                                            for k, v in value.items():
                                                                if type(v) is dict:
                                                                    file_w.write(str(k) + ':')
                                                                    file_w.write('\n')
                                                                    # Change key to value and value to key
                                                                    for ki, vi in v.items():
                                                                        if vi in cell_types_dict:
                                                                            cell_types_dict[vi] += ki + '; '
                                                                        else:
                                                                            cell_types_dict[vi] = ki + '; '
                                                                    for ki, vi in cell_types_dict.items():
                                                                        file_w.write('\t' + str(ki) + ' = ' + str(vi))
                                                                        file_w.write('\n')

                                                except HTTPError as http_error:
                                                    print(f'HTTP error occurred: {http_error}')
                                                except Exception as error:
                                                    print(f'Other error occurred: {error}')
                                                file_w.write('\n')
                                                file_w.write('\n')
                                                sleep(2)

                                    # -------------------------------------------
                                    # Ensembl Cross Reference
                                    # -------------------------------------------

                                    # Looking for id first gene
                                    if choice_Cross_References == 'Y':
                                        attributes = non_overlap_list[l1][8]
                                        if 'ID=gene' in attributes:
                                            start_position = attributes.index(":") + 1
                                            if start_position > 20:
                                                start_position = attributes.index("-") + 1

                                            end_position = attributes.find(';', start_position)
                                            id_search = attributes[start_position:end_position]

                                            file_w.write(
                                                '------------------------------------------------------------------------------------')
                                            file_w.write('\n')
                                            file_w.write('Cross References : ' + id_search)
                                            file_w.write('\n')
                                            file_w.write(
                                                '------------------------------------------------------------------------------------')

                                            file_w.write('\n')
                                            # Mounting query URL in the Ensembl XReference
                                            server = "https://rest.ensembl.org"
                                            ext = f'/xrefs/id/{id_search}?all_levels=1'
                                            r = requests.get(server + ext, headers={"Content-Type": "application/json"})
                                            print()
                                            if not r.ok:
                                                r.raise_for_status()
                                            decoded = r.json()
                                            for value in decoded:
                                                if type(value) is dict:
                                                    for k, v in value.items():
                                                        if type(v) is list and len(v) == 0:
                                                            pass
                                                        elif v == '':
                                                            pass
                                                        else:
                                                            file_w.write("'" + str(k) + "'" + ":" + "'" + str(v) + "'" + ', ')
                                                    file_w.write('.')
                                                    file_w.write('\n')
                                            file_w.write('\n')
                                            sleep(1)
                                            file_w.write('\n')
                                            file_w.write('\n')

                                        # Looking for id second gene
                                        attributes = ''
                                        attributes = non_overlap_list[l2][8]
                                        if 'ID=gene' in attributes:
                                            start_position = attributes.index(":") + 1
                                            if start_position > 20:
                                                start_position = attributes.index("-") + 1
                                            end_position = attributes.find(';', start_position)
                                            id_search = attributes[start_position:end_position]

                                            file_w.write(
                                                '------------------------------------------------------------------------------------')
                                            file_w.write('\n')
                                            file_w.write('Cross References : ' + id_search)
                                            file_w.write('\n')
                                            file_w.write(
                                                '------------------------------------------------------------------------------------')
                                            file_w.write('\n')

                                            server = "https://rest.ensembl.org"
                                            ext = f'/xrefs/id/{id_search}?all_levels=1'
                                            r = requests.get(server + ext, headers={"Content-Type": "application/json"})
                                            if not r.ok:
                                                r.raise_for_status()
                                            decoded = r.json()
                                            for value in decoded:
                                                if type(value) is dict:
                                                    for k, v in value.items():
                                                        if type(v) is list and len(v) == 0:
                                                            pass
                                                        elif v == '':
                                                            pass
                                                        else:
                                                            file_w.write(f"'{str(k)}': '{str(v)}', ")
                                                    file_w.write('.')
                                                    file_w.write('\n')
                                            file_w.write('\n')
                                            file_w.write('\n')
                                    sleep(1)
                                    cont_interval += 1
                        except Exception as err:
                            print(err)


    except Exception as erro:
        print(erro)
    finally:
        file_w.write('\n')
        file_w.write('Number of intergenic intervals: ' + str(cont_interval))
        file_w.write('\n')
        file_w.close()

    print('\n')
    print('--' * 15)
    print('RESULT')
    print('--' * 15)
    print(f'File generated with the results: {name_file_write}')
    print('=-' * 15)


# ---------------------------------------------------------------------------------
#                                   Program Start
# ---------------------------------------------------------------------------------

def main():
    option = ''
    plt.ioff()
    while True:
        file_r, file_name, interval = define_files()

        list_gene = parse_gff(file_r)

        non_overlap_list, total_amount_overlap = remove_overlap(list_gene)

        ordered_dictionary = assemble_table(non_overlap_list, interval, len(list_gene), total_amount_overlap)

        while True:
            plot = str(input('Deseja plotar a tabela de Genes ? (Y/ N):')).strip().upper()[0]
            if plot in 'YN':
                break
            else:
                print('Error, answer only Y or N.')

        if plot == 'Y':
            interval_plot = 6
            interval_plot = input('Informe o intervalo da plotagem (default 6) :')

            plot_genes(ordered_dictionary, int(interval_plot), file_name)

        # option_fasta = None
        # while True:
        #     option_fasta = str(input('Do you want to get the FASTA sequence? (Y/N) ')).strip().upper()[0]
        #     if option_fasta in 'YN':
        #         break
        #     else:
        #         print('Error, answer only Y or N.')
        # if option_fasta == 'Y':

        fasta(non_overlap_list, file_name)

        while True:
            option = str(input('Do you want to analyze another GFF file? (Y/N) ')).strip().upper()[0]
            if option in 'YN':
                break
            else:
                print('Error, answer only Y or N.')
        if option == 'N':
            print('End of Program! Have a Good One.')
            break
        else:
            input(os.name)
            if os.name == "nt":
                os.system("cls")
            else:
                os.system("clear")


main()
