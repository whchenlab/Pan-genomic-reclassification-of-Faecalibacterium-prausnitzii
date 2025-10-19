# -*- coding: iso-8859-1 -*-

import csv
import os

def count_occurrences(row):
    """Í³¼Æ»ùÒòÔÚ¶àÉÙ¸ö¾úÖêÖÐ³öÏÖ"""
    count = 0
    # ´ÓµÚ15ÁÐ¿ªÊ¼Í³¼Æ£¨Ë÷Òý14£¬ÒòÎªPython´Ó0¿ªÊ¼£©
    for value in row[14:]:
        if value.strip() != '':
            count += 1
    return count

def find_first_occurrence(row, strains):
    """ÕÒµ½»ùÒòµÚÒ»¸ö³öÏÖµÄ¾úÖê¼°Æä×ùÎ»ºÅ"""
    # ´ÓµÚ15ÁÐ¿ªÊ¼²éÕÒ£¨Ë÷Òý14£©
    for i, value in enumerate(row[14:], start=14):
        if value.strip() != '':
            strain_index = i - 14  # ¾úÖêÁÐ±íµÄË÷Òý
            return strains[strain_index], value.strip()
    return None, None

def extract_gene_sequence(faa_file, locus_tag):
    """´ÓFAAÎÄ¼þÖÐÌáÈ¡Ö¸¶¨×ùÎ»ºÅµÄ»ùÒòÐòÁÐ"""
    if not os.path.exists(faa_file):
        print(f"¾¯¸æ: ÎÄ¼þ {faa_file} ²»´æÔÚ")
        return None
    
    sequence = []
    capture = False
    
    with open(faa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # ¼ì²éÊÇ·ñÊÇÄ¿±ê»ùÒòµÄ×¢ÊÍÐÐ
                gene_id = line.split()[0][1:]  # ÌáÈ¡>ºóÃæµ½µÚÒ»¸ö¿Õ¸ñÇ°µÄÄÚÈÝ
                if gene_id == locus_tag:
                    capture = True
                    sequence.append(line.strip())
                else:
                    capture = False
            elif capture:
                sequence.append(line.strip())
    
    if not sequence:
        print(f"¾¯¸æ: ÔÚ {faa_file} ÖÐÎ´ÕÒµ½×ùÎ»ºÅ {locus_tag} µÄ»ùÒò")
        return None
    
    return '\n'.join(sequence)

def main():
    # ÊäÈëÎÄ¼þÂ·¾¶
    csv_file = './roaryresult/gene_presence_absence.csv'
    faa_dir = './faas/'
    
    # ¼ì²éÊäÈëÎÄ¼þÊÇ·ñ´æÔÚ
    if not os.path.exists(csv_file):
        print(f"´íÎó: ÎÄ¼þ {csv_file} ²»´æÔÚ")
        return
    
    if not os.path.exists(faa_dir):
        print(f"´íÎó: Ä¿Â¼ {faa_dir} ²»´æÔÚ")
        return
    
    # ´ò¿ªÊä³öÎÄ¼þ
    core_file = open('core.faa', 'w')
    shell_file = open('shell.faa', 'w')
    cloud_file = open('cloud.faa', 'w')
    
    try:
        with open(csv_file, 'r') as f:
            reader = csv.reader(f)
            # ¶ÁÈ¡±íÍ·
            header = next(reader)
            
            # 1. ¼ÆËã×Ü¾úÖêÊý£¨´ÓµÚ15ÁÐ¿ªÊ¼£©
            total_strains = len(header) - 14  # 15ÁÐ¶ÔÓ¦Ë÷Òý14
            print(f"×Ü¾úÖêÊý: {total_strains}")
            
            # ¼ÆËã¸÷ÀàÐÍ»ùÒòµÄãÐÖµ
            core_threshold = total_strains * 0.95
            shell_low_threshold = total_strains * 0.15
            
            print(f"Core»ùÒòãÐÖµ: ¡Ý {core_threshold:.1f} ¸ö¾úÖê")
            print(f"Shell»ùÒòãÐÖµ: {shell_low_threshold:.1f} - {core_threshold:.1f} ¸ö¾úÖê")
            print(f"Cloud»ùÒòãÐÖµ: < {shell_low_threshold:.1f} ¸ö¾úÖê")
            
            # »ñÈ¡¾úÖêÃû³ÆÁÐ±í£¨´ÓµÚ15ÁÐ¿ªÊ¼£©
            strains = header[14:]
            
            # ´¦ÀíÃ¿¸ö»ùÒò
            gene_count = 0
            for row in reader:
                gene_count += 1
                if gene_count % 100 == 0:
                    print(f"ÒÑ´¦Àí {gene_count} ¸ö»ùÒò...")
                
                # »ùÒòÃû³ÆÔÚµÚÒ»ÁÐ
                gene_name = row[0]
                
                # Í³¼Æ¸Ã»ùÒòÔÚ¶àÉÙ¾úÖêÖÐ³öÏÖ
                occurrence_count = count_occurrences(row)
                
                # È·¶¨»ùÒòÀàÐÍ
                if occurrence_count >= core_threshold:
                    gene_type = 'core'
                    output_file = core_file
                elif occurrence_count >= shell_low_threshold:
                    gene_type = 'shell'
                    output_file = shell_file
                else:
                    gene_type = 'cloud'
                    output_file = cloud_file
                
                # ÕÒµ½µÚÒ»¸ö³öÏÖ¸Ã»ùÒòµÄ¾úÖê¼°Æä×ùÎ»ºÅ
                strain_name, locus_tag = find_first_occurrence(row, strains)
                if not strain_name or not locus_tag:
                    print(f"¾¯¸æ: »ùÒò {gene_name} ÔÚËùÓÐ¾úÖêÖÐ¾ùÎ´³öÏÖ£¬Ìø¹ý")
                    continue
                
                # ¹¹½¨FAAÎÄ¼þÂ·¾¶
                faa_file = os.path.join(faa_dir, f"{strain_name}.faa")
                
                # ÌáÈ¡»ùÒòÐòÁÐ
                sequence = extract_gene_sequence(faa_file, locus_tag)
                if sequence:
                    # Ð´Èë¶ÔÓ¦µÄÊä³öÎÄ¼þ
                    output_file.write(f"{sequence}\n")
        
        print(f"´¦ÀíÍê³É£¬¹²´¦Àí {gene_count} ¸ö»ùÒò")
        print("½á¹ûÒÑ·Ö±ðÐ´Èë core.faa, shell.faa ºÍ cloud.faa")
    
    finally:
        # È·±£ÎÄ¼þ±»¹Ø±Õ
        core_file.close()
        shell_file.close()
        cloud_file.close()

if __name__ == "__main__":
    main()
