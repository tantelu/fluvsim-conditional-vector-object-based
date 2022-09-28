#include "filehelper.h";
#include "channel.h";
#include "omp.h"
#include "mathhelp.h"
#include <set>

void movecenterline() {
	vector<corssinfo> infos;
	int wellid = 5;
	auto pt = readwell("\\data\\singlechcondition\\well\\well.txt");
	auto wells = pt.get();
	if (wells != nullptr) {
		double dis = 1;
		double wellspace = 50;
		channelzone zone(0, Coordinate(-10, 50, 0), 550);
		modelpar modelpar(1, 0.25, 0.5, -99);
		channelpar par(zone, modelpar, 25, 5, 100, 50);
		channel ch(par);

		vector<int> c_ig;
		vector<facie> c_fs;
		c_fs.reserve(wells->get_size());
		for (size_t i = 0; i < wells->get_size(); i++)
		{
			c_fs.emplace_back((facie)wells->get_at(i).get_facie());
		}
		wellsfilter filter(c_ig, c_fs);
		vector<double> ratios;

		ch.build_center_nodes(nullptr);
		{
			auto centerpt = ch.get_centerline();
			auto center = centerpt.get();
			writepoints("\\data\\singlechcondition\\line\\base.txt", center);
		}
		{
			auto sectionpt = ch.get_cross_chsectionline(&wells->get_with_id(wellid).get_pos(), dis);
			auto section = sectionpt.get();
			if (section != nullptr) {
				string outpath = "\\data\\singlechcondition\\section\\base.txt";
				writepoints(outpath, section);
			}
		}
		{
			corssinfo info;
			ch.well_cross_info(*wells, wellspace / 2.0, info.alllen, info.cross, info.sand, info.mud);
			infos.push_back(info);
		}
		for (size_t i = 1; i <= 5; i++)
		{
			int stepNums = 5;

			ch.move_center_nodes(wells, &wellspace, &dis, &stepNums);

			auto centerpt = ch.get_centerline();
			auto center = centerpt.get();
			if (center != nullptr) {
				string outpath = "\\data\\singlechcondition\\line\\";
				outpath.append(std::to_string(i));
				outpath.append(".txt");
				writepoints(outpath, center);
			}
			auto sectionpt = ch.get_cross_chsectionline(&wells->get_with_id(wellid).get_pos(), dis);
			auto section = sectionpt.get();
			if (section != nullptr) {
				string outpath = "\\data\\singlechcondition\\section\\";
				outpath.append(std::to_string(i));
				outpath.append(".txt");
				writepoints(outpath, section);
			}

			corssinfo info;
			ch.well_cross_info(*wells, wellspace / 2.0, info.alllen, info.cross, info.sand, info.mud);
			infos.push_back(info);
		}

		writeinfos("\\data\\singlechcondition\\ratio.txt", infos);
	}
}

void selectchswithdifferentlen() {

	vector<int> sandids;
	auto pt = readwell("\\data\\singlechcondition\\well\\well.txt");
	auto wells = pt.get();
	for (size_t i = 0; i < wells->get_size(); i++)
	{
		if ((facie)wells->get_at(i).get_facie() == (facie)1) {
			sandids.push_back(wells->get_at(i).get_id());
		}
	}
	int groupnums = 10;

	string compareoutfloder = "\\data\\chlencompare\\channels\\";
	clearfolder(compareoutfloder);
	for (int len = 100; len <= 500; len += 50)
	{
		
		string lenoutfolder = compareoutfloder;
		lenoutfolder.append(std::to_string(len));
		createfolder(lenoutfolder);

		string channelfolder = "\\data\\channels\\";
		channelfolder.append(std::to_string(len));
		vector<string> allfiles;
		find(channelfolder, allfiles);
		unordered_map<int, vector<string>> idmaps;
		classify_file_by_id(allfiles, idmaps);
		for (size_t groupi = 0; groupi < groupnums; groupi++)
		{
			string groupifolder = lenoutfolder;
			groupifolder.append("\\");
			groupifolder.append(std::to_string(groupi));
			createfolder(groupifolder);

			set<string> needfiles2;
			for (auto it = idmaps.begin(); it != idmaps.end(); it++)
			{
				set<string> temp;
				auto list = (*it).second;
				int minselect = min((int)(list.size() - 1), 2);
				while (temp.size() < minselect) {
					int shift = (int)geo_math::random(0, list.size() - 0.01);
					temp.insert(list[shift]);
				}
				for (auto it2 = temp.begin(); it2 != temp.end(); it2++)
				{
					needfiles2.insert(*it2);
				}
			}
			vector<string> needfiles;
			needfiles.assign(needfiles2.begin(), needfiles2.end());

			vector<vector<int>> ch_has_sandids;
			vector<int> ch_index_in_allchs;
			set<int> all;
			ch_has_sandids.reserve(needfiles.size());
			for (size_t i = 0; i < needfiles.size(); i++)
			{
				vector<string> res;
				boost::split(res, needfiles[i], boost::is_any_of(",."), boost::token_compress_on);
				if (res.size() > 2) {
					vector<int> ids;
					for (size_t j = 1; j < res.size() - 1; j++)
					{
						int welli = std::stoi(res[j]);
						ids.push_back(welli);
						all.insert(welli);
					}
					ch_has_sandids.push_back(ids);
					ch_index_in_allchs.push_back(i);
				}
			}
			vector<int> all2;
			all2.assign(all.begin(), all.end());
			auto res = geo_math::smallestSufficientTeam(all2, ch_has_sandids);
			
			vector<string> selects;
			selects.reserve(res.size());
			for (size_t i = 0; i < res.size(); i++)
			{
				string chfile = channelfolder;
				chfile.append("\\");
				chfile.append(needfiles[ch_index_in_allchs[res[i]]]);
				selects.push_back(chfile);
			}
			for (size_t i = 0; i < selects.size(); i++)
			{
				int pos = selects[i].find_last_of('\\');
				string s(selects[i].substr(pos + 1));
				string target = groupifolder;
				target.append("\\");
				target.append(s);
				copyfile(selects[i], target);
			}
		}
	}

}


void chs_to_model() {
	Coordinate gridsize(1, 1, 0.25);
	modelpar modelpar(gridsize.x, gridsize.z, gridsize.x / 2, -99);

	int groupnums = 10;
	string compareoutfloder = "\\data\\chlencompare\\channels\\";
	string modelfloder = "\\data\\chlencompare\\models\\";
	for (int len = 100; len <= 500; len += 50)
	{
		string lenoutfolder = compareoutfloder;
		lenoutfolder.append(std::to_string(len));
		for (size_t groupi = 0; groupi < groupnums; groupi++)
		{
			grid2d grid(Envelope(0, 500, 0, 300), gridsize);
			string groupifolder = lenoutfolder;
			groupifolder.append("\\");
			groupifolder.append(std::to_string(groupi));
			modelcenterstomodels(groupifolder, modelpar, grid);

			string gslibfile = modelfloder;
			gslibfile.append(std::to_string(len));
			gslibfile.append("-");
			gslibfile.append(std::to_string(groupi));
			gslibfile.append("-");
			gslibfile.append("model.txt");
			writegslib(gslibfile, grid);
		}
	}
}

void generate_ch_database() {
	rect3d rect3(0, 500, 0, 300, -10, 0);
	Coordinate gridsize(1, 1, 0.25);
	modelpar modelpar(gridsize.x, gridsize.z, gridsize.x / 2, -99);

	auto pt = readwell("\\data\\channels\\well\\well.txt");
	auto wells = pt.get();

	vector<well3d> sandwells;
	for (size_t welli = 0; welli < wells->get_size(); welli++)
	{
		if (wells->get_at(welli).get_facie() == facie(1)) {
			sandwells.push_back(wells->get_at(welli));
		}
	}

	size_t maxcount = 500;
	int everywellmincount = 10;
	double dis = 1;
	double wellspace = 50;
	int sandStemNums = 10;
	int checkindex = min(0.9 * maxcount, maxcount - 1.0 * everywellmincount * sandwells.size());//开始检查井点条件化的起始索引
	for (int len = 500; len <= 500; len += 50)
	{
		vector<int> ch_add_count(sandwells.size(), 0);
		string folder = "\\data\\channels\\";
		folder.append(std::to_string(len));
		clearfolder(folder);
		createfolder(folder);
		for (int index = 0; index < maxcount; index++)
		{
			double startdotx = geo_math::random(-10, 450.0 - len);
			double startdoty = geo_math::random(rect3.get_ymin(), rect3.get_ymax());
			double amp = geo_math::random(20, 35);
			double wave = geo_math::random(50, 100);
			double angel = geo_math::random(-0.1, 0.1);
			double width = geo_math::random(20, 35);
			double thick = geo_math::random(4, 5);
			channelzone zone(angel, Coordinate(startdotx, startdoty, 0), min(1.0 * len, 550.0 - startdotx));
			channelpar par(zone, modelpar, width, thick, wave, amp);
			channel ch(par);
			ch.build_center_nodes(nullptr);
			ch.move_center_nodes(wells, &wellspace, &dis, &sandStemNums);
			string file = "\\data\\channels\\";
			file.append(std::to_string(len));
			file.append("\\");
			file.append(std::to_string(index));
			for (size_t welli = 0; welli < sandwells.size(); welli++)
			{
				if (ch.is_well_match(sandwells[welli])) {
					ch_add_count[welli]++;
					file.append(",");
					file.append(std::to_string(sandwells[welli].get_id()));
				}
			}
			file.append(".txt");
			if (file.find(',') != string::npos) {
				ch.save_to_file(file);
			}
			else {
				index--;
			}
		}
		for (size_t i = 0; i < sandwells.size(); i++)
		{
			int chindex = maxcount;
			while (ch_add_count[i] < everywellmincount) {
				double width = geo_math::random(20, 35);
				double amp = geo_math::random(0, width * 2);
				double wave = geo_math::random(50, 100);
				double angel = 0;
				double thick = geo_math::random(4, 5);
				double startdotx = startdotx = sandwells[i].get_pos().x - geo_math::random(len * 0.25, len * 0.75);
				double startdoty = startdoty = sandwells[i].get_pos().y + geo_math::random(-width, width);
				channelzone zone(angel, Coordinate(startdotx, startdoty, 0), len);
				channelpar par(zone, modelpar, width, thick, wave, amp);
				channel ch(par);
				ch.build_center_nodes(nullptr);
				ch.move_center_nodes(wells, &wellspace, &dis, &sandStemNums);
				if (ch.is_well_match(sandwells[i])) {
					ch_add_count[i]++;
					string file = "\\data\\channels\\";
					file.append(std::to_string(len));
					file.append("\\");
					file.append(std::to_string(chindex));
					file.append(",");
					file.append(std::to_string(sandwells[i].get_id()));
					file.append(".txt");
					ch.save_to_file(file);
					chindex++;
				}
			}
		}
	}
}

void nonstationaritydatabase() {
	auto pt = readwell("\\data\\nonstationarity\\wells\\well.txt");
	auto wells = pt.get();
	double wellspace = 50;
	double dis = 1;
	int sandStemNums = 15;
	vector<well3d> sandwells;
	for (size_t welli = 0; welli < wells->get_size(); welli++)
	{
		if (wells->get_at(welli).get_facie() == facie(1)) {
			sandwells.push_back(wells->get_at(welli));
		}
	}

	Coordinate gridsize(1, 1, 0.25);
	Envelope range(0, 400, 0, 400);
	LineSegment top(Coordinate(0, 400), Coordinate(400, 400));
	LineSegment right(Coordinate(400, 0), Coordinate(400, 400));
	modelpar modelpar(1, 0.25, 0.5, -99);
	string path = "\\data\\nonstationarity\\centers\\";
	clearfolder(path);
	createfolder(path);
	vector<int> ch_add_count(sandwells.size(), 0);
	for (size_t i = 0; i < 500; i++)
	{
		double len = 400 * 1.414;
		double width = 30;
		double thick = 6;
		double angel = geo_math::random(0, 1.57);
		Coordinate start = i % 2 == 0 ? Coordinate(geo_math::random(0, 300), 0) : Coordinate(0, geo_math::random(0, 300));
		width = width - start.distance(Coordinate(0, 0)) / 100 * 0.1;
		Coordinate end(start.x + cos(angel) * 400, start.y + sin(angel) * 400);
		LineSegment line(start, end);
		auto topcross = line.intersection(top);
		auto rightcross = line.intersection(right);
		if (!topcross.isNull()) {
			len = topcross.distance(start);
		}
		if (!rightcross.isNull()) {
			len = rightcross.distance(start);
		}
		channelzone subzone(angel, start, len);
		channelpar par(subzone, modelpar, width, thick, 40, 20);
		channel subch(par);
		subch.build_center_nodes(0.1);
		subch.move_center_nodes(wells, &wellspace, &dis, &sandStemNums);

		string linefile = path;
		linefile.append(std::to_string(i + 1));

		for (size_t welli = 0; welli < sandwells.size(); welli++)
		{
			if (subch.is_well_match(sandwells[welli])) {
				ch_add_count[welli]++;
				linefile.append(",");
				linefile.append(std::to_string(sandwells[welli].get_id()));
			}
		}
		linefile.append(".txt");
		if (linefile.find(',') != string::npos) {
			subch.save_to_file(linefile);
		}
		else {
			i--;
		}
	}
	/*grid2d grid(Envelope(0, 500, 0, 500), gridsize);
	writegridvalue(grid, chs);
	string gslib = "\\data\\nonstationarity\\models\\model.txt";
	writegslib(gslib, grid);*/
}

void nonstationarityselect() {
	modelpar modelpar(1, 0.25, 0.5, -99);
	auto pt = readwell("\\data\\nonstationarity\\wells\\well.txt");
	auto wells = pt.get();
	vector<well3d> sandwells;
	for (size_t welli = 0; welli < wells->get_size(); welli++)
	{
		if (wells->get_at(welli).get_facie() == facie(1)) {
			sandwells.push_back(wells->get_at(welli));
		}
	}
	Coordinate orginal(0, 0);
	std::sort(sandwells.begin(), sandwells.end(), [orginal](const well3d& well1, const well3d& well2)->bool {return well1.get_pos().distanceSquared(orginal) < well2.get_pos().distanceSquared(orginal); });
	set<int> hasids;

	string floder = "\\data\\nonstationarity\\centers\\";
	string selectfolder = "\\data\\nonstationarity\\selectfolder\\";
	clearfolder(selectfolder);
	createfolder(selectfolder);
	vector<string> allfiles;
	find(floder, allfiles);
	random_shuffle(allfiles.begin(), allfiles.end());
	vector<channel> chs;
	for (size_t i = 0; i < sandwells.size(); i++)
	{
		auto id = sandwells[i].get_id();
		if (hasids.find(id) == hasids.end()) {
			bool has = false;
			for (size_t j = 0; j < allfiles.size(); j++)
			{
				vector<string> res;
				boost::split(res, allfiles[j], boost::is_any_of(",."), boost::token_compress_on);
				if (res.size() > 2) {
					for (size_t k = 1; k < res.size() - 1; k++)
					{
						int welli = std::stoi(res[k]);
						if (welli == id) {
							has = true;
						}
					}
				}
				if (has) {
					string file = floder;
					file.append(allfiles[j]);
					channel ch(file, modelpar);
					ch.cut_center_with_others(chs, sandwells[i]);
					for (size_t k = 0; k < sandwells.size(); k++)
					{
						if (ch.is_well_match(sandwells[k])) {
							hasids.insert(sandwells[k].get_id());
						}
					}
					chs.push_back(ch);
					string des = selectfolder;
					des.append(allfiles[j]);
					ch.save_to_file(des);
					break;
				}
			}
		}
	}
	Coordinate gridsize(1, 1, 0.25);
	grid2d grid(Envelope(0, 400, 0, 400), gridsize);
	writegridvalue(grid, chs);
	string gslib = "\\data\\nonstationarity\\models\\model.txt";
	writegslib(gslib, grid);
}

void stationarityselectwithcut() {
	srand(4545);
	modelpar modelpar(1, 0.25, 0.5, -99);
	auto pt = readwell("\\data\\channels\\well\\well.txt");
	auto wells = pt.get();
	vector<well3d> sandwells;
	for (size_t welli = 0; welli < wells->get_size(); welli++)
	{
		if (wells->get_at(welli).get_facie() == facie(1)) {
			sandwells.push_back(wells->get_at(welli));
		}
	}
	Coordinate orginal(0, 0);
	std::sort(sandwells.begin(), sandwells.end(), [orginal](const well3d& well1, const well3d& well2)->bool {return well1.get_pos().distanceSquared(orginal) < well2.get_pos().distanceSquared(orginal); });
	set<int> hasids;

	string floder = "\\data\\channels\\500\\";
	string cutselectfolder = "\\data\\cutcomparenocut\\cut\\";
	clearfolder(cutselectfolder);
	createfolder(cutselectfolder);
	vector<string> allfiles;
	find(floder, allfiles);
	random_shuffle(allfiles.begin(), allfiles.end());
	vector<channel> cutchs;
	for (size_t i = 0; i < sandwells.size(); i++)
	{
		auto id = sandwells[i].get_id();
		if (hasids.find(id) == hasids.end()) {
			bool has = false;
			for (size_t j = 0; j < allfiles.size(); j++)
			{
				vector<string> res;
				boost::split(res, allfiles[j], boost::is_any_of(",."), boost::token_compress_on);
				if (res.size() > 2) {
					for (size_t k = 1; k < res.size() - 1; k++)
					{
						int welli = std::stoi(res[k]);
						if (welli == id) {
							has = true;
						}
					}
				}
				if (has) {
					string file = floder;
					file.append(allfiles[j]);
					channel ch(file, modelpar);

					ch.cut_center_with_others(cutchs, sandwells[i]);
					for (size_t k = 0; k < sandwells.size(); k++)
					{
						if (ch.is_well_match(sandwells[k])) {
							hasids.insert(sandwells[k].get_id());
						}
					}
					cutchs.push_back(ch);
					string cutdes = cutselectfolder;
					cutdes.append(allfiles[j]);
					ch.save_to_file(cutdes);
					break;
				}
			}
		}
	}
	/*Coordinate gridsize(1, 1, 0.25);
	grid2d grid(Envelope(0, 400, 0, 400), gridsize);
	writegridvalue(grid, chs);
	string gslib = "\\data\\nonstationarity\\models\\model.txt";
	writegslib(gslib, grid);*/
}

void savechtomodel(string folder, string des) {
	Coordinate gridsize(1, 1, 0.25);
	modelpar modelpar(gridsize.x, gridsize.z, gridsize.x / 2, -99);
	grid2d grid(Envelope(0, 500, 0, 300), gridsize);
	modelcenterstomodels(folder, modelpar, grid);
	writegslib(des, grid);
}

void genereal_nostationary() {
	Coordinate gridsize(1, 1, 0.25);
	modelpar modelpar(gridsize.x, gridsize.z, gridsize.x / 2, -99);
	Envelope env = Envelope(0, 250, 0, 250);
	grid2d grid(env, gridsize);
	string folder = "\\data\\nonstationarity";
	vector<channel> chs;
	for (size_t i = 0; i < 1; i++)
	{
		string gslibfile = folder.append(std::to_string(i));
		for (size_t j = 0; j <= 10; j++)
		{
			double angel = -1+j*0.2;
			double width = geo_math::random(6, 9);
			double wave = geo_math::random(50, 60);
			double amp = geo_math::random(10,15);
			channelzone zone(angel, Coordinate(-5, 125, 0), 800);
			channelpar par(zone, modelpar, width, 5, wave, amp);
			channel ch(par);
			ch.build_center_nodes(nullptr);
			ch.get_in_channel_pos(&grid);
			chs.push_back(ch);
		}
		writegslib(gslibfile, grid);
	}
}

int main()
{
	genereal_nostationary();
	//chs_to_model();
	//generate_ch_database();
	//nonstationaritydatabase();
	//stationarityselectwithcut();
	/*string folder = "\\data\\cutcomparenocut\\cut\\";
	string des = "\\data\\cutcomparenocut\\cut.txt";
	savechtomodel(folder, des);*/

	std::cout << "Hello World!\n";
}


