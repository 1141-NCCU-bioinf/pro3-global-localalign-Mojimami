import numpy as np
import pandas as pd

def alignment(input_path, score_path, output_path, aln, gap):
    headers = []
    sequences = []
    with open(input_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                headers.append(line)
            elif line:
                sequences.append(line)
    
    seq_x = sequences[0]
    seq_y = sequences[1]
    
    matrix_df = pd.read_csv(score_path, sep='\s+', comment='#', index_col=0)
    score_map = matrix_df.to_dict()
    
    n_rows = len(seq_y) + 1
    n_cols = len(seq_x) + 1
    
    scores = np.zeros((n_rows, n_cols))
    paths = [[set() for _ in range(n_cols)] for _ in range(n_rows)]
    
    if aln == 'global':
        scores[0, :] = np.arange(n_cols) * gap
        scores[:, 0] = np.arange(n_rows) * gap
        
        for r in range(1, n_rows):
            paths[r][0].add('v')
        for c in range(1, n_cols):
            paths[0][c].add('h')
        
        for r in range(1, n_rows):
            for c in range(1, n_cols):
                diag = scores[r-1, c-1] + score_map[seq_x[c-1]][seq_y[r-1]]
                vert = scores[r-1, c] + gap
                horz = scores[r, c-1] + gap
                
                best = max(diag, vert, horz)
                scores[r, c] = best
                
                if diag == best:
                    paths[r][c].add('d')
                if vert == best:
                    paths[r][c].add('v')
                if horz == best:
                    paths[r][c].add('h')
        
        r, c = n_rows - 1, n_cols - 1
        align_x = ''
        align_y = ''
        
        while r > 0 or c > 0:
            if r == 0:
                align_x = seq_x[c-1] + align_x
                align_y = '-' + align_y
                c -= 1
            elif c == 0:
                align_x = '-' + align_x
                align_y = seq_y[r-1] + align_y
                r -= 1
            else:
                val = scores[r, c]
                if val == scores[r-1, c-1] + score_map[seq_x[c-1]][seq_y[r-1]]:
                    align_x = seq_x[c-1] + align_x
                    align_y = seq_y[r-1] + align_y
                    r -= 1
                    c -= 1
                elif val == scores[r-1, c] + gap:
                    align_x = '-' + align_x
                    align_y = seq_y[r-1] + align_y
                    r -= 1
                else:
                    align_x = seq_x[c-1] + align_x
                    align_y = '-' + align_y
                    c -= 1
        
        with open(output_path, 'w') as f:
            f.write(headers[0] + '\n')
            f.write(align_x + '\n')
            f.write(headers[1] + '\n')
            f.write(align_y + '\n')
    
    else:
        peak = 0
        peak_cells = []
        
        for r in range(1, n_rows):
            for c in range(1, n_cols):
                diag = scores[r-1, c-1] + score_map[seq_x[c-1]][seq_y[r-1]]
                vert = scores[r-1, c] + gap
                horz = scores[r, c-1] + gap
                zero = 0
                
                options = {'d': diag, 'v': vert, 'h': horz, 'z': zero}
                best = max(options.values())
                scores[r, c] = best
                
                for key, val in options.items():
                    if val == best:
                        paths[r][c].add(key)
                
                if best > peak:
                    peak = best
                    peak_cells = [(r, c)]
                elif best == peak and best > 0:
                    peak_cells.append((r, c))
        
        found = []
        
        def backtrack(row, col, ax, ay):
            if row <= 0 and col <= 0:
                found.append((ax, ay))
                return
            elif col <= 0:
                backtrack(row-1, col, '-' + ax, seq_y[row-1] + ay)
                return
            elif row <= 0:
                backtrack(row, col-1, seq_x[col-1] + ax, '-' + ay)
                return
            
            if 'z' in paths[row][col] or scores[row, col] == 0:
                ax = seq_x[col-1] + ax
                ay = seq_y[row-1] + ay
                found.append((ax, ay))
                return
            
            if 'd' in paths[row][col]:
                backtrack(row-1, col-1, seq_x[col-1] + ax, seq_y[row-1] + ay)
            if 'h' in paths[row][col]:
                backtrack(row, col-1, seq_x[col-1] + ax, '-' + ay)
            if 'v' in paths[row][col]:
                backtrack(row-1, col, '-' + ax, seq_y[row-1] + ay)
        
        for start_r, start_c in peak_cells:
            backtrack(start_r, start_c, '', '')
        
        if found:
            longest = max(len(p[0]) for p in found)
            found = [(px, py) for px, py in found if len(px) == longest]
            found = sorted(set(found), key=lambda pair: (pair[0], pair[1]))
            
            with open(output_path, 'w') as f:
                for align_x, align_y in found:
                    f.write(headers[0] + '\n')
                    f.write(align_x + '\n')
                    f.write(headers[1] + '\n')
                    f.write(align_y + '\n')
        else:
            open(output_path, 'w').close()

if __name__ == "__main__":
    alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/result_global.fasta", "global", -10)
    alignment("examples/test_local.fasta", "examples/pam100.txt", "examples/result_local.fasta", "local", -10)
