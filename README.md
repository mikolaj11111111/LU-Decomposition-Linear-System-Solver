# LU Decomposition Linear System Solver

## Description
Solves systems of linear equations using LU decomposition with partial pivoting. Factors matrix A into lower triangular L and upper triangular U matrices, then solves using forward and backward substitution.

## Features
- **LU Decomposition** - Factorizes A = PLU with partial pivoting
- **Partial Pivoting** - Improves numerical stability through row swapping
- **Forward Substitution** - Solves Lz = Pb for intermediate vector z
- **Backward Substitution** - Solves Ux = z for final solution x
- **Solution Verification** - Checks accuracy by computing A*x vs b
- **Step-by-Step Visualization** - Shows decomposition process iteration by iteration

## Mathematical Background
- **LU Decomposition**: A = PLU where P is permutation, L is lower triangular, U is upper triangular
- **Partial Pivoting**: Selects largest pivot to minimize numerical errors
- **Two-Stage Solution**: Lz = Pb → Ux = z gives solution to Ax = b

## Input Format
Files should contain:
```
N = [matrix size]
b: [vector values]
A: [matrix rows]
```

## Usage
Processes multiple test files automatically: `LU_gr1INO 1.txt`, `gauss_elimination_gr1IO_A.txt`, etc.
Shows complete decomposition process with intermediate steps and final verification.
