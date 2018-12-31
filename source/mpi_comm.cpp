#include "mpi_comm.hpp"
#include <array>
#ifdef RICH_MPI

namespace
{
	std::array<double, 7> Primitive2Array(Primitive const& cell)
	{
		std::array<double, 7> res;
		res[0] = cell.density;
		res[1] = cell.energy;
		res[2] = cell.entropy;
		res[3] = cell.LastCool;
		res[4] = cell.velocity;
		res[5] = cell.pressure;
		res[6] = static_cast<double>(cell.sticker);
		return res;
	}

	Primitive Array2Primitive(std::array<double, 7> const& data)
	{
		Primitive res;
		res.density = data[0];
		res.energy = data[1];
		res.entropy = data[2];
		res.LastCool = data[3];
		res.velocity = data[4];
		res.pressure = data[5];
		res.sticker = static_cast<unsigned char>(data[6]);
		return res;
	}

	std::array<double, 5> Extensive2Array(Extensive const& cell)
	{
		std::array<double, 5> res;
		res[0] = cell.mass;
		res[1] = cell.energy;
		res[2] = cell.entropy;
		res[3] = cell.momentum;
		res[4] = cell.et;
		return res;
	}

	Extensive Array2Extensive(std::array<double, 5> const& data)
	{
		Extensive res;
		res.mass = data[0];
		res.energy = data[1];
		res.entropy = data[2];
		res.momentum = data[3];
		res.et = data[4];
		return res;
	}

	std::vector<double> ExtensiveVec2Double(std::vector<Extensive> const& cells)
	{
		size_t N = cells.size();
		std::vector<double> res;
		res.reserve(5 * N);
		for (size_t i = 0; i < N; ++i)
		{
			std::array<double, 5> temp = Extensive2Array(cells[i]);
			for (size_t j = 0; j < 5; ++j)
				res.push_back(temp[j]);
		}
		return res;
	}

	std::vector<Extensive> VecDouble2Extensive(std::vector<double> const& data)
	{
		size_t N = data.size() / 5;
		std::vector<Extensive> res(N);
		for (size_t i = 0; i < N; ++i)
		{
			std::array<double, 5> temp;
			for (size_t j = 0; j < 5; ++j)
				temp[j] = data[i * 5 + j];
			res[i] = Array2Extensive(temp);
		}
		return res;
	}

	std::vector<double> PrimitiveVec2Double(std::vector<Primitive> const& cells)
	{
		size_t N = cells.size();
		std::vector<double> res;
		res.reserve(7 * N);
		for (size_t i = 0; i < N; ++i)
		{
			std::array<double, 7> temp = Primitive2Array(cells[i]);
			for (size_t j = 0; j < 7; ++j)
				res.push_back(temp[j]);
		}
		return res;
	}

	std::vector<Primitive> VecDouble2Primitive(std::vector<double> const& data)
	{
		size_t N = data.size() / 7;
		std::vector<Primitive> res(N);
		for (size_t i = 0; i < N; ++i)
		{
			std::array<double, 7> temp;
			for (size_t j = 0; j < 7; ++j)
				temp[j] = data[i * 7 + j];
			res[i] = Array2Primitive(temp);
		}
		return res;
	}
}

std::array<Primitive, 4> SendRecvPrimitive(std::vector<Primitive> const& cells)
{
	int rank = 0,ws=0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	std::array<double, 14> data,data2,recv;
	std::array<double, 7> temp;
	std::array<Primitive, 4> res;
	MPI_Request req;
	if (rank == 0)
	{
		size_t N = cells.size() - 1;
		temp = Primitive2Array(cells[N-1]);
		std::copy(temp.begin(), temp.end(), data.begin());
		temp = Primitive2Array(cells[N]);
		std::copy(temp.begin(), temp.end(), data.begin()+7);
		MPI_Isend(&data[0], 14, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD,&req);
		MPI_Recv(&recv[0], 14, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		std::copy(recv.begin(), recv.begin() + 7, temp.begin());
		res[2] = Array2Primitive(temp);
		std::copy(recv.begin() + 7, recv.end() , temp.begin());
		res[3] = Array2Primitive(temp);
	}
	else
	{
		if (rank == (ws - 1))
		{
			temp = Primitive2Array(cells[0]);
			std::copy(temp.begin(), temp.end(), data.begin());
			temp = Primitive2Array(cells[1]);
			std::copy(temp.begin(), temp.end(), data.begin() + 7);
			MPI_Isend(&data[0], 14, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &req);
			MPI_Recv(&recv[0], 14, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			std::copy(recv.begin(), recv.begin() + 7, temp.begin());
			res[0] = Array2Primitive(temp);
			std::copy(recv.begin() + 7, recv.end(), temp.begin());
			res[1] = Array2Primitive(temp);
		}
		else
		{
			MPI_Status status;
			size_t N = cells.size() - 1;
			temp = Primitive2Array(cells[N - 1]);
			std::copy(temp.begin(), temp.end(), data2.begin());
			temp = Primitive2Array(cells[N]);
			std::copy(temp.begin(), temp.end(), data2.begin() + 7);
			temp = Primitive2Array(cells[0]);
			std::copy(temp.begin(), temp.end(), data.begin());
			temp = Primitive2Array(cells[1]);
			std::copy(temp.begin(), temp.end(), data.begin() + 7);
			MPI_Isend(&data[0], 14, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &req);
			MPI_Isend(&data2[0], 14, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &req);
			MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&recv[0], 14, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (status.MPI_SOURCE == (rank - 1))
			{
				std::copy(recv.begin(), recv.begin() + 7, temp.begin());
				res[0] = Array2Primitive(temp);
				std::copy(recv.begin() + 7, recv.end(), temp.begin());
				res[1] = Array2Primitive(temp);
			}
			else
			{
				std::copy(recv.begin(), recv.begin() + 7, temp.begin());
				res[2] = Array2Primitive(temp);
				std::copy(recv.begin() + 7, recv.end(), temp.begin());
				res[3] = Array2Primitive(temp);
			}
			MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&recv[0], 14, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (status.MPI_SOURCE == (rank - 1))
			{
				std::copy(recv.begin(), recv.begin() + 7, temp.begin());
				res[0] = Array2Primitive(temp);
				std::copy(recv.begin() + 7, recv.end(), temp.begin());
				res[1] = Array2Primitive(temp);
			}
			else
			{
				std::copy(recv.begin(), recv.begin() + 7, temp.begin());
				res[2] = Array2Primitive(temp);
				std::copy(recv.begin() + 7, recv.end(), temp.begin());
				res[3] = Array2Primitive(temp);
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return res;
}

std::array<double,4> SendRecvEdges(std::vector<double> const & edges)
{
	int rank = 0, ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	std::array<double, 4> res;
	std::array<double, 2> recv;
	MPI_Request req;
	if (rank == 0)
	{
		size_t N = edges.size() - 2;
		MPI_Isend(&edges[N-1], 2, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &req);
		MPI_Recv(&recv[0], 2, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		res[2] = recv[0];
		res[3] = recv[1];
	}
	else
	{
		if (rank == (ws - 1))
		{
			MPI_Isend(&edges[1], 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &req);
			MPI_Recv(&recv[0],2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			res[0] = recv[0];
			res[1] = recv[1];
		}
		else
		{
			MPI_Status status;
			size_t N = edges.size() - 2;
			MPI_Isend(&edges[1], 2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &req);
			MPI_Isend(&edges[N-1], 2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &req);
			MPI_Probe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&recv[0], 2, MPI_DOUBLE, status.MPI_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (status.MPI_SOURCE == (rank - 1))
			{
				res[0] = recv[0];
				res[1] = recv[1];
			}
			else
			{
				res[2] = recv[0];
				res[3] = recv[1];
			}
			MPI_Probe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&recv[0], 2, MPI_DOUBLE, status.MPI_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (status.MPI_SOURCE == (rank - 1))
			{
				res[0] = recv[0];
				res[1] = recv[1];
			}
			else
			{
				res[2] = recv[0];
				res[3] = recv[1];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return res;
}

void RedistributeExtensives(std::vector<Extensive> &cells, std::vector<double> &edges)
{
	int nlocal = cells.size();
	int ntotal = 0;
	MPI_Allreduce(&nlocal, &ntotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	int rank = 0, ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	size_t index_lower = (rank * ntotal) / ws;
	size_t index_upper = ((rank+1) * ntotal) / ws -1;
	std::vector<int> nperproc(ws, 0),disp(ws,0),newn(ws,0);
	nlocal *= 5;
	MPI_Gather(&nlocal, 1, MPI_INT, &nperproc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	std::vector<double> tosend = ExtensiveVec2Double(cells);
	std::vector<double> torecv(2);
	if (rank == 0)
	{
		torecv.resize(ntotal * 5);
		for (size_t i = 1; i < ws; ++i)
			disp[i] = nperproc[i - 1] + disp[i - 1];
	}
	MPI_Gatherv(&tosend[0], tosend.size(), MPI_DOUBLE, &torecv[0], &nperproc[0], &disp[0], MPI_DOUBLE, 0,
		MPI_COMM_WORLD);
	tosend.resize(5 * (index_upper - index_lower + 1));
	// send the data
	if (rank == 0)
	{
		for (size_t i = 0; i < ws; ++i)
			newn[i] = 5*(((i + 1) * ntotal) / ws - (i * ntotal) / ws);
		for (size_t i = 1; i < ws; ++i)
			disp[i] = newn[i - 1] + disp[i - 1];
	}
	MPI_Scatterv(&torecv[0], &newn[0], &disp[0], MPI_DOUBLE, &tosend[0],
		5*(index_upper - index_lower + 1),MPI_DOUBLE, 0, MPI_COMM_WORLD);
	cells = VecDouble2Extensive(tosend);
	// deal with the edges
	if (rank == 0)
	{
		torecv.resize(ntotal + 1);
		for (size_t i = 0; i < ws; ++i)
			nperproc[i] /= 5;
		for (size_t i = 1; i < ws; ++i)
			disp[i] = nperproc[i - 1] + disp[i - 1];
	}
	MPI_Gatherv(&edges[1], edges.size() - 1, MPI_DOUBLE, &torecv[1], &nperproc[0], &disp[0],
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(rank==0)
		torecv[0] = edges[0];
	tosend.resize(index_upper - index_lower + 2);
	if (rank == 0)
	{
		for (size_t i = 0; i < ws; ++i)
		{
			newn[i] /= 5;
			newn[i]++;
			if(i>0)
				disp[i] =newn[i-1] + disp[i-1] - 1;
		}
	}
	MPI_Scatterv(&torecv[0], &newn[0], &disp[0], MPI_DOUBLE, &tosend[0],
		(index_upper - index_lower + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	edges = tosend;
}
void ConsolidateData(std::vector<Primitive>& cells, std::vector<double>& edges, std::vector<std::vector<double> >& append)
{
	int nlocal = cells.size();
	int ntotal = 0;
	MPI_Allreduce(&nlocal, &ntotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	int rank = 0, ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	std::vector<int> nperproc(ws, 0), disp(ws, 0);
	nlocal *= 7;
	MPI_Gather(&nlocal, 1, MPI_INT, &nperproc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	std::vector<double> tosend = PrimitiveVec2Double(cells);
	std::vector<double> torecv(2);
	if (rank == 0)
	{
		torecv.resize(ntotal * 7);
		for (size_t i = 1; i < ws; ++i)
			disp[i] = nperproc[i - 1] + disp[i - 1];
	}
	MPI_Gatherv(&tosend[0], tosend.size(), MPI_DOUBLE, &torecv[0], &nperproc[0], &disp[0], MPI_DOUBLE, 0,
		MPI_COMM_WORLD);
	if (rank == 0)
		cells = VecDouble2Primitive(torecv);
	// deal with the edges
	if (rank == 0)
	{
		torecv.resize(ntotal + 1);
		for (size_t i = 0; i < nperproc.size(); ++i)
			nperproc[i] /= 7;
		for (size_t i = 1; i < ws; ++i)
			disp[i] = nperproc[i - 1] + disp[i - 1];
	}
	MPI_Gatherv(&edges[1], edges.size() - 1, MPI_DOUBLE, &torecv[1], &nperproc[0], &disp[0],
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		torecv[0] = edges[0];
		edges = torecv;
	}
	// deal with appendices
	for (size_t i = 0; i < append.size(); ++i)
	{
		if (rank == 0)
		{
			torecv.resize(ntotal);
		}
		MPI_Gatherv(&append[i][0], append[i].size(), MPI_DOUBLE, &torecv[0], &nperproc[0], &disp[0],
			MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (rank == 0)
			append[i] = torecv;
	}
}
#endif