import numpy as np
import pandas as pd

def alignment(input_path, score_path, output_path, aln, gap):
    """Main alignment function"""
    
    # Read FASTA file
    with open(input_path, 'r') as f:
        lines = f.readlines()
    headers = [lines[i].strip() for i in range(0, len(lines), 2)]
    seqs = [lines[i].strip() for i in range(1, len(lines), 2)]
    s1, s2 = seqs[0], seqs[1]
    
    # Load scoring matrix using pandas
    df = pd.read_csv(score_path, sep='\s+', comment='#', index_col=0)
    score_dict = df.to_dict()
    
    # Initialize matrices
    m, n = len(s2) + 1, len(s1) + 1
    dp = np.zeros((m, n))
    
    if aln == 'global':
        # Global alignment initialization
        dp[0, :] = np.arange(n) * gap
        dp[:, 0] = np.arange(m) * gap
        
        # Fill matrix
        for i in range(1, m):
            for j in range(1, n):
                match = dp[i-1, j-1] + score_dict.get(s1[j-1], {}).get(s2[i-1], gap)
                delete = dp[i-1, j] + gap
                insert = dp[i, j-1] + gap
                dp[i, j] = max(match, delete, insert)
        
        # Traceback from bottom-right
        i, j = m-1, n-1
        a1, a2 = '', ''
        
        while i > 0 or j > 0:
            if i == 0:
                a1 = s1[j-1] + a1
                a2 = '-' + a2
                j -= 1
            elif j == 0:
                a1 = '-' + a1
                a2 = s2[i-1] + a2
                i -= 1
            else:
                current = dp[i, j]
                if current == dp[i-1, j-1] + score_dict.get(s1[j-1], {}).get(s2[i-1], gap):
                    a1 = s1[j-1] + a1
                    a2 = s2[i-1] + a2
                    i -= 1
                    j -= 1
                elif current == dp[i-1, j] + gap:
                    a1 = '-' + a1
                    a2 = s2[i-1] + a2
                    i -= 1
                else:
                    a1 = s1[j-1] + a1
                    a2 = '-' + a2
                    j -= 1
        
        results = [(a1, a2)]
        
    else:  # local alignment
        # Keep track of paths at each cell
        paths = [[[] for _ in range(n)] for _ in range(m)]
        
        # Fill matrix and track maximum
        max_score = 0
        max_pos = []
        
        for i in range(1, m):
            for j in range(1, n):
                match = dp[i-1, j-1] + score_dict.get(s1[j-1], {}).get(s2[i-1], gap)
                delete = dp[i-1, j] + gap
                insert = dp[i, j-1] + gap
                
                score = max(0, match, delete, insert)
                dp[i, j] = score
                
                if score > 0:
                    if score == match:
                        paths[i][j].append('m')
                    if score == delete:
                        paths[i][j].append('d')
                    if score == insert:
                        paths[i][j].append('i')
                
                if score > max_score:
                    max_score = score
                    max_pos = [(i, j)]
                elif score == max_score and score > 0:
                    max_pos.append((i, j))
        
        # Traceback all paths
        all_alns = []
        
        def trace(i, j, a1='', a2=''):
            if dp[i, j] == 0:
                if a1:
                    all_alns.append((a1, a2))
                return
            
            for move in paths[i][j]:
                if move == 'm':
                    trace(i-1, j-1, s1[j-1] + a1, s2[i-1] + a2)
                elif move == 'd':
                    trace(i-1, j, '-' + a1, s2[i-1] + a2)
                elif move == 'i':
                    trace(i, j-1, s1[j-1] + a1, '-' + a2)
        
        for pos in max_pos:
            trace(pos[0], pos[1])
        
        # Filter by length and sort
        if all_alns:
            max_len = max(len(x[0]) for x in all_alns)
            results = [x for x in all_alns if len(x[0]) == max_len]
            results = sorted(set(results))
        else:
            results = []
    
    # Write output
    with open(output_path, 'w') as f:
        for a1, a2 in results:
            f.write(f'{headers[0]}\n{a1}\n{headers[1]}\n{a2}\n')

# Test the implementation
if __name__ == "__main__":
    alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/result_global.fasta", "global", -10)
    alignment("examples/test_local.fasta", "examples/pam100.txt", "examples/result_local.fasta", "local", -10)
