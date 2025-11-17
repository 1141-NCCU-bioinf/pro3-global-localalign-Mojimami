import numpy as np
import pandas as pd

def alignment(input_path, score_path, output_path, aln, gap):
    sequences = []
    headers = []
    with open(input_path, 'r') as file:
        current_seq = ""
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                headers.append(line)
                current_seq = ""
            else:
                current_seq += line
        if current_seq:
            sequences.append(current_seq)
    
    seq_a = sequences[0]
    seq_b = sequences[1]
    
    scoring_matrix = pd.read_csv(score_path, sep='\s+', comment='#', index_col=0)
    
    score_lookup = {}
    for row_aa in scoring_matrix.index:
        score_lookup[row_aa] = {}
        for col_aa in scoring_matrix.columns:
            score_lookup[row_aa][col_aa] = int(scoring_matrix.loc[row_aa, col_aa])
    
    rows = len(seq_b) + 1
    cols = len(seq_a) + 1
    
    score_matrix = np.zeros((rows, cols), dtype=float)
    
    direction_matrix = [[[] for _ in range(cols)] for _ in range(rows)]
    
    if aln == 'global':
        for col_idx in range(cols):
            score_matrix[0][col_idx] = gap * col_idx
            if col_idx > 0:
                direction_matrix[0][col_idx] = ['left']
        
        for row_idx in range(rows):
            score_matrix[row_idx][0] = gap * row_idx
            if row_idx > 0:
                direction_matrix[row_idx][0] = ['up']
        
        for row_idx in range(1, rows):
            for col_idx in range(1, cols):
                char_a = seq_a[col_idx - 1]
                char_b = seq_b[row_idx - 1]
                
                diagonal_score = score_matrix[row_idx-1][col_idx-1] + score_lookup[char_b][char_a]
                up_score = score_matrix[row_idx-1][col_idx] + gap
                left_score = score_matrix[row_idx][col_idx-1] + gap
                
                max_score = max(diagonal_score, up_score, left_score)
                score_matrix[row_idx][col_idx] = max_score
                
                if diagonal_score == max_score:
                    direction_matrix[row_idx][col_idx].append('diag')
                if up_score == max_score:
                    direction_matrix[row_idx][col_idx].append('up')
                if left_score == max_score:
                    direction_matrix[row_idx][col_idx].append('left')
        
        aligned_a = ""
        aligned_b = ""
        row_pos = rows - 1
        col_pos = cols - 1
        
        while row_pos > 0 or col_pos > 0:
            if row_pos == 0:
                aligned_a = seq_a[col_pos - 1] + aligned_a
                aligned_b = '-' + aligned_b
                col_pos -= 1
            elif col_pos == 0:
                aligned_a = '-' + aligned_a
                aligned_b = seq_b[row_pos - 1] + aligned_b
                row_pos -= 1
            else:
                current_score = score_matrix[row_pos][col_pos]
                char_a = seq_a[col_pos - 1]
                char_b = seq_b[row_pos - 1]
                
                if current_score == score_matrix[row_pos-1][col_pos-1] + score_lookup[char_b][char_a]:
                    aligned_a = char_a + aligned_a
                    aligned_b = char_b + aligned_b
                    row_pos -= 1
                    col_pos -= 1
                elif current_score == score_matrix[row_pos-1][col_pos] + gap:
                    aligned_a = '-' + aligned_a
                    aligned_b = char_b + aligned_b
                    row_pos -= 1
                else:
                    aligned_a = char_a + aligned_a
                    aligned_b = '-' + aligned_b
                    col_pos -= 1
        
        with open(output_path, 'w') as out_file:
            out_file.write(headers[0] + '\n')
            out_file.write(aligned_a + '\n')
            out_file.write(headers[1] + '\n')
            out_file.write(aligned_b + '\n')
    
    else:
        max_score_value = 0
        max_score_positions = []
        
        for row_idx in range(1, rows):
            for col_idx in range(1, cols):
                char_a = seq_a[col_idx - 1]
                char_b = seq_b[row_idx - 1]
                
                diagonal_score = score_matrix[row_idx-1][col_idx-1] + score_lookup[char_b][char_a]
                up_score = score_matrix[row_idx-1][col_idx] + gap
                left_score = score_matrix[row_idx][col_idx-1] + gap
                
                max_score = max(0, diagonal_score, up_score, left_score)
                score_matrix[row_idx][col_idx] = max_score
                
                if max_score > 0:
                    if diagonal_score == max_score:
                        direction_matrix[row_idx][col_idx].append('diag')
                    if up_score == max_score:
                        direction_matrix[row_idx][col_idx].append('up')
                    if left_score == max_score:
                        direction_matrix[row_idx][col_idx].append('left')
                elif max_score == 0:
                    if diagonal_score == 0:
                        direction_matrix[row_idx][col_idx].append('diag')
                    if up_score == 0:
                        direction_matrix[row_idx][col_idx].append('up')
                    if left_score == 0:
                        direction_matrix[row_idx][col_idx].append('left')
                
                if max_score > max_score_value:
                    max_score_value = max_score
                    max_score_positions = [(row_idx, col_idx)]
                elif max_score == max_score_value and max_score > 0:
                    max_score_positions.append((row_idx, col_idx))
        
        all_alignments = set()
        
        def traceback_recursive(r_pos, c_pos, aln_a, aln_b):
            if r_pos == 0 or c_pos == 0:
                if aln_a:
                    all_alignments.add((aln_a, aln_b))
                return
            
            current_score = score_matrix[r_pos][c_pos]
            
            if current_score == 0:
                for direction in direction_matrix[r_pos][c_pos]:
                    if direction == 'diag' and r_pos > 0 and c_pos > 0:
                        prev_score = score_matrix[r_pos-1][c_pos-1]
                        if prev_score > 0:
                            new_a = seq_a[c_pos - 1] + aln_a
                            new_b = seq_b[r_pos - 1] + aln_b
                            all_alignments.add((new_a, new_b))
                    elif direction == 'up' and r_pos > 0:
                        prev_score = score_matrix[r_pos-1][c_pos]
                        if prev_score > 0:
                            new_a = '-' + aln_a
                            new_b = seq_b[r_pos - 1] + aln_b
                            all_alignments.add((new_a, new_b))
                    elif direction == 'left' and c_pos > 0:
                        prev_score = score_matrix[r_pos][c_pos-1]
                        if prev_score > 0:
                            new_a = seq_a[c_pos - 1] + aln_a
                            new_b = '-' + aln_b
                            all_alignments.add((new_a, new_b))
                if not direction_matrix[r_pos][c_pos] and aln_a:
                    all_alignments.add((aln_a, aln_b))
                return
            
            for direction in direction_matrix[r_pos][c_pos]:
                if direction == 'diag':
                    new_a = seq_a[c_pos - 1] + aln_a
                    new_b = seq_b[r_pos - 1] + aln_b
                    traceback_recursive(r_pos - 1, c_pos - 1, new_a, new_b)
                elif direction == 'up':
                    new_a = '-' + aln_a
                    new_b = seq_b[r_pos - 1] + aln_b
                    traceback_recursive(r_pos - 1, c_pos, new_a, new_b)
                elif direction == 'left':
                    new_a = seq_a[c_pos - 1] + aln_a
                    new_b = '-' + aln_b
                    traceback_recursive(r_pos, c_pos - 1, new_a, new_b)
        
        for start_r, start_c in max_score_positions:
            traceback_recursive(start_r, start_c, "", "")
        
        alignment_list = list(all_alignments)
        if alignment_list:
            max_length = max(len(aln[0]) for aln in alignment_list)
            filtered_alignments = [(a, b) for a, b in alignment_list if len(a) == max_length]
            
            filtered_alignments.sort(key=lambda x: (x[0], x[1]))
            
            with open(output_path, 'w') as out_file:
                for aligned_a, aligned_b in filtered_alignments:
                    out_file.write(headers[0] + '\n')
                    out_file.write(aligned_a + '\n')
                    out_file.write(headers[1] + '\n')
                    out_file.write(aligned_b + '\n')
        else:
            with open(output_path, 'w') as out_file:
                pass

if __name__ == "__main__":
    alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/result_global.fasta", "global", -10)
    alignment("examples/test_local.fasta", "examples/pam100.txt", "examples/result_local.fasta", "local", -10)
