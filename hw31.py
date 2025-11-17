import numpy as np
import pandas as pd

def alignment(input_path, score_path, output_path, aln, gap):
    # Parse FASTA file
    with open(input_path, 'r') as file:
        lines = file.readlines()
    
    headers, sequences = [], []
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('>'):
            headers.append(line)
        else:
            sequences.append(line)
    
    seq1, seq2 = sequences[0], sequences[1]
    header1, header2 = headers[0], headers[1]
    
    # Load substitution matrix
    score_matrix = pd.read_csv(score_path, sep='\s+', comment='#', index_col=0)
    
    # Create lookup dictionary
    substitution_scores = {}
    for residue1 in score_matrix.index:
        substitution_scores[residue1] = {}
        for residue2 in score_matrix.columns:
            substitution_scores[residue1][residue2] = int(score_matrix.at[residue1, residue2])
    
    # Initialize DP matrices
    m, n = len(seq2) + 1, len(seq1) + 1
    dp_score = np.zeros((m, n), dtype=int)
    dp_path = [[[] for _ in range(n)] for _ in range(m)]
    
    if aln == 'global':
        # Initialize first row and column
        for j in range(n):
            dp_score[0][j] = gap * j
            if j > 0:
                dp_path[0][j].append('left')
        
        for i in range(m):
            dp_score[i][0] = gap * i
            if i > 0:
                dp_path[i][0].append('up')
        
        # Fill DP table
        for i in range(1, m):
            for j in range(1, n):
                base1 = seq1[j-1]
                base2 = seq2[i-1]
                
                score_diag = dp_score[i-1][j-1] + substitution_scores[base2][base1]
                score_up = dp_score[i-1][j] + gap
                score_left = dp_score[i][j-1] + gap
                
                max_score = max(score_diag, score_up, score_left)
                dp_score[i][j] = max_score
                
                # Store all paths that achieve max score
                if score_diag == max_score:
                    dp_path[i][j].append('diag')
                if score_up == max_score:
                    dp_path[i][j].append('up')
                if score_left == max_score:
                    dp_path[i][j].append('left')
        
        # Traceback for global alignment
        aligned_seq1 = ""
        aligned_seq2 = ""
        i, j = m - 1, n - 1
        
        while i > 0 or j > 0:
            if i == 0:
                aligned_seq1 = seq1[j-1] + aligned_seq1
                aligned_seq2 = '-' + aligned_seq2
                j -= 1
            elif j == 0:
                aligned_seq1 = '-' + aligned_seq1
                aligned_seq2 = seq2[i-1] + aligned_seq2
                i -= 1
            else:
                current_score = dp_score[i][j]
                base1 = seq1[j-1]
                base2 = seq2[i-1]
                
                # Choose path based on score calculation
                if current_score == dp_score[i-1][j-1] + substitution_scores[base2][base1]:
                    aligned_seq1 = base1 + aligned_seq1
                    aligned_seq2 = base2 + aligned_seq2
                    i, j = i - 1, j - 1
                elif current_score == dp_score[i-1][j] + gap:
                    aligned_seq1 = '-' + aligned_seq1
                    aligned_seq2 = base2 + aligned_seq2
                    i -= 1
                else:
                    aligned_seq1 = base1 + aligned_seq1
                    aligned_seq2 = '-' + aligned_seq2
                    j -= 1
        
        # Write result
        with open(output_path, 'w') as out:
            out.write(header1 + '\n')
            out.write(aligned_seq1 + '\n')
            out.write(header2 + '\n')
            out.write(aligned_seq2 + '\n')
    
    else:  # local alignment
        # Track maximum score positions
        max_score_value = 0
        max_positions = []
        
        # Fill DP table with modified logic
        for i in range(1, m):
            for j in range(1, n):
                base1 = seq1[j-1]
                base2 = seq2[i-1]
                
                # Calculate all possible scores
                match_score = dp_score[i-1][j-1] + substitution_scores[base2][base1]
                delete_score = dp_score[i-1][j] + gap
                insert_score = dp_score[i][j-1] + gap
                
                # Local alignment allows 0 as minimum
                best_score = max(0, match_score, delete_score, insert_score)
                dp_score[i][j] = best_score
                
                # Track which moves gave the best score
                moves = []
                if best_score == match_score and best_score > 0:
                    moves.append('diag')
                if best_score == delete_score and best_score > 0:
                    moves.append('up')
                if best_score == insert_score and best_score > 0:
                    moves.append('left')
                if best_score == 0:
                    if match_score == 0:
                        moves.append('diag')
                    if delete_score == 0:
                        moves.append('up')
                    if insert_score == 0:
                        moves.append('left')
                
                dp_path[i][j] = moves
                
                # Update maximum positions
                if best_score > max_score_value:
                    max_score_value = best_score
                    max_positions = [(i, j)]
                elif best_score == max_score_value and best_score > 0:
                    max_positions.append((i, j))
        
        # Collect all alignments using recursive traceback
        alignments = []
        
        def traceback(row, col, align1, align2):
            # Check boundaries first
            if row <= 0 and col <= 0:
                if align1:
                    alignments.append((align1, align2))
                return
            
            if row <= 0:
                traceback(row, col-1, seq1[col-1] + align1, '-' + align2)
                return
            
            if col <= 0:
                traceback(row-1, col, '-' + align1, seq2[row-1] + align2)
                return
            
            # Check if we've reached the end of alignment
            if dp_score[row][col] == 0:
                # Include the boundary characters
                final_align1 = seq1[col-1] + align1
                final_align2 = seq2[row-1] + align2
                alignments.append((final_align1, final_align2))
                return
            
            # Continue traceback through all possible paths
            for move in dp_path[row][col]:
                if move == 'diag':
                    traceback(row-1, col-1, seq1[col-1] + align1, seq2[row-1] + align2)
                elif move == 'up':
                    traceback(row-1, col, '-' + align1, seq2[row-1] + align2)
                elif move == 'left':
                    traceback(row, col-1, seq1[col-1] + align1, '-' + align2)
        
        # Start traceback from all maximum positions
        for pos_i, pos_j in max_positions:
            traceback(pos_i, pos_j, '', '')
        
        # Process results
        if alignments:
            # Remove duplicates
            unique_alignments = list(set(alignments))
            
            # Filter by maximum length
            max_length = max(len(a[0]) for a in unique_alignments)
            filtered = [(a1, a2) for a1, a2 in unique_alignments if len(a1) == max_length]
            
            # Sort lexicographically
            filtered.sort(key=lambda x: (x[0], x[1]))
            
            # Write all alignments
            with open(output_path, 'w') as out:
                for aligned1, aligned2 in filtered:
                    out.write(header1 + '\n')
                    out.write(aligned1 + '\n')
                    out.write(header2 + '\n')
                    out.write(aligned2 + '\n')
        else:
            # Create empty file if no alignments found
            open(output_path, 'w').close()

if __name__ == "__main__":
    alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/result_global.fasta", "global", -10)
    alignment("examples/test_local.fasta", "examples/pam100.txt", "examples/result_local.fasta", "local", -10)