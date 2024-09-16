/**
 * @Author: Your name
 * @Date:   2024-09-16 11:02:48
 * @Last Modified by:   Your name
 * @Last Modified time: 2024-09-16 15:01:00
 */
#include"head.h"
#include"LBSolver.h"
int main()
{
	int Nx = 64;
	int Ny = 64;
	int Nz = 64;
	LBSolver LB(Nx, Ny, Nz);
	LB.set_para();
	LB.solver();
	return 0;
}