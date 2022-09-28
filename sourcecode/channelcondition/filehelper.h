#pragma once
#include <iostream>
#include <fstream>
#include "channel.h";
#include <stdio.h>
#include <Windows.h>
#include<unordered_map>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

struct corssinfo
{
	double alllen = 0;
	double cross = 0;
	double sand = 0;
	double mud = 0;
};

unique_ptr<well3d_sequence> readwell(string path) {
	ifstream infile;
	infile.open(path, ios::in);
	if (!infile.is_open())
	{
		return unique_ptr<well3d_sequence>(nullptr);
	}
	vector<well3d> wells;
	string buf;
	int id = 0;
	while (getline(infile, buf))
	{
		vector<string> res;
		boost::split(res, buf, boost::is_any_of(",."), boost::token_compress_on);
		if (res.size() == 4) {
			auto x = std::stod(res[0]);
			auto y = std::stod(res[1]);
			auto thick = std::stod(res[3]);
			auto f = std::stoi(res[2]);
			auto bot = 0.0 - thick;
			wells.push_back(well3d(id, Coordinate(x, y), (facie)f, 0.0, bot));
		}
		id++;
	}
	return unique_ptr<well3d_sequence>(new well3d_sequence(&wells));
}

int writepoints(const string& path, const CoordinateArraySequence* center) {
	ofstream out_file(path);

	if (out_file.is_open()) {
		if (center != nullptr) {
			for (size_t i = 0; i < center->getSize(); i++)
			{
				auto cor = center->getAt(i);
				out_file << cor.x << ',' << cor.y << endl;
			}
		}
	}
	else {
		return 0;
	}
	out_file.close();
}

int writeinfos(const string& path, const vector<corssinfo>& vs) {
	ofstream out_file(path);

	if (out_file.is_open()) {
		for (size_t i = 0; i < vs.size(); i++)
		{
			corssinfo v = vs[i];
			out_file << v.alllen << " " << v.cross << " " << v.sand << " " << v.mud << endl;
		}
	}
	else {
		return 0;
	}
	out_file.close();
}

void writegslib(const string& path, const grid3d& grid) {
	ofstream out_file(path);
	if (out_file.is_open()) {
		auto size = grid.get_size();
		auto range = grid.get_range();
		auto xcount = grid.get_xcount();
		auto ycount = grid.get_ycount();
		auto zcount = grid.get_zcount();
		out_file << "petrel(" << grid.get_xcount() << 'x' << grid.get_ycount() << 'x' << grid.get_zcount() << ')' << endl;
		out_file << "1" << endl;
		out_file << "fluvsim" << endl;
		auto vs = grid.get_facies();
		for (int i = 0; i < grid.get_gridunits_count(); i++)
		{
			out_file << *(vs + i) << endl;
		}
	}
	out_file.close();
}

void writegridvalue(grid3d& grid, const vector<channel>& channels) {
	for (size_t i = 0; i < channels.size(); i++)
	{
		channels.at(i).get_in_channel_pos(&grid);
	}
}

void writegslib(const string& path, const grid2d& grid) {
	ofstream out_file(path);
	if (out_file.is_open()) {
		/*auto size = grid.get_size();
		auto range = grid.get_range();
		auto xcount = grid.get_xcount();
		auto ycount = grid.get_ycount();
		out_file << "petrel(" << grid.get_xcount() << 'x' << grid.get_ycount() << ')' << endl;
		out_file << "1" << endl;
		out_file << "fluvsim" << endl;*/
		auto vs = grid.get_facies();
		for (int i = 0; i < grid.get_gridunits_count(); i++)
		{
			out_file << vs[i] << endl;
		}
	}
	out_file.close();
}

void writegridvalue(grid2d& grid, const vector<channel>& channels) {
	for (size_t i = 0; i < channels.size(); i++)
	{
		channels.at(i).get_in_channel_pos(&grid);
	}
}

void find(string lpPath, vector<string>& fileList)
{
	char szFind[MAX_PATH];
	WIN32_FIND_DATA FindFileData;

	strcpy(szFind, lpPath.c_str());
	strcat(szFind, "\\*.txt");

	HANDLE hFind = ::FindFirstFile(szFind, &FindFileData);
	if (INVALID_HANDLE_VALUE == hFind)    return;

	while (true)
	{
		if (FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
		{
			if (FindFileData.cFileName[0] != '.')
			{
				char szFile[MAX_PATH];
				strcpy(szFile, lpPath.c_str());
				strcat(szFile, "\\");
				strcat(szFile, (char*)(FindFileData.cFileName));
				find(szFile, fileList);
			}
		}
		else
		{
			//std::cout << FindFileData.cFileName << std::endl;
			fileList.push_back(FindFileData.cFileName);
		}
		if (!FindNextFile(hFind, &FindFileData))    break;
	}
	FindClose(hFind);
}

int copyfile(const string& source, const string& target) {
	ifstream in(source, ios::binary);
	ofstream out(target, ios::binary);
	if (!in) {
		printf("open file error");
		return -1;
	}
	if (!out) {
		printf("open file error");
		return -1;
	}
	out << in.rdbuf();
	in.close();
	out.close();
	return 1;
}\

void classify_file_by_id(const vector<string>& filenames, unordered_map<int, vector<string>>& idmaps) {
	for (auto it = filenames.begin(); it != filenames.end(); it++)
	{
		vector<string> res;
		boost::split(res, (*it), boost::is_any_of(",."), boost::token_compress_on);
		if (res.size() > 2) {
			vector<int> ids;
			for (size_t j = 1; j < res.size() - 1; j++)
			{
				int welli = std::stoi(res[j]);
				auto map = idmaps.find(welli);
				if (map != idmaps.end()) {
					map->second.push_back(*it);
				}
				else
				{
					vector<string> tempfile;
					tempfile.push_back(*it);
					idmaps.insert(pair<int, vector<string>>(welli, tempfile));
				}
			}
		}
	}
}

void clearfolder(const string& floder) {
	string flodercmd = "rd /s /q ";
	flodercmd.append(floder);
	system(flodercmd.c_str());
}

void createfolder(const string& floder) {
	string command;
	command = "mkdir " + floder;
	system(command.c_str());
}

void modelcenterstomodels(const string& floder, const modelpar& modelpar,grid2d& grid) {
	auto pos = floder.find_last_of('\\');
	vector<string> allfiles;
	find(floder, allfiles);
	for (size_t i = 0; i < allfiles.size(); i++)
	{
		string tempfile = floder;
		if (pos != floder.length() - 1) {
			tempfile.append("\\");
		}
		tempfile.append(allfiles[i]);
		channel ch(tempfile, modelpar);
		ch.get_in_channel_pos(&grid);
	}
}


