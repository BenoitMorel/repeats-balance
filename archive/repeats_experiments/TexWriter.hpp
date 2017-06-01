#ifndef _TEXWRITER_H_
#define _TEXWRITER_H_

#include <fstream>
#include <vector>
#include <string>

template<typename T>
struct Plot {
  struct HorizontalLine {
    T y;
    bool dashed;
    std::string color;
    std::string caption;
  };

  struct Legend {
    std::string color;
    std::string text;
  };

  struct Arrow {
    double x;
    T ymin;
    T ymax;
    std::string color;
  };

  Plot(T max, const std::string caption) {
    this->max = max;
    this->caption = caption;
    this->histo = false;
  }

  
  Plot &add_plot(const std::vector<T> &plot, const std::string &color = "blue", char mark = '.') {
    plots.push_back(&plot);
    plots_colors.push_back(color);
    plots_marks.push_back(mark);
    return *this;
  }

  Plot &is_histo(bool histogram) {
    histo = histogram;
    return *this;
  }

  Plot &xlabel(const std::string &label) {
    xcaption = label;
    return *this;
  }

  Plot &ylabel(const std::string &label) {
    ycaption = label;
    return *this;
  }

  Plot &add_arrow(double x, T ymin, T ymax, const std::string &color = "red") {
    Arrow a;
    a.x = x;
    a.ymin = ymin;
    a.ymax = ymax;
    a.color = color;
    arrows.push_back(a);
    return *this;
  }

  Plot &add_legend(const std::string &color, const std::string &text) {
    Legend l;
    l.color = color;
    l.text = text;
    legend.push_back(l);
    return *this;
  }

  Plot &add_horizontal_line(T y, const std::string &color = "red", 
                            bool dashed=true, const std::string &caption = "") {
    HorizontalLine line;
    line.y = y;
    line.color = color;
    line.dashed = dashed;
    line.caption = caption;
    lines.push_back(line);
    return *this;
  }

  friend std::ostream& operator<< (std::ostream &os, const Plot<T> &p) {
    if (p.plots.size() == 0) {
      std::cout << "Error, no plot to plot !" << std::endl;
    }
    os << "\\begin{minipage}{0.49\\textwidth}" << std::endl;
    os << "\\begin{tikzpicture}[scale=0.75]" << std::endl;
    if (p.histo) {
      os << "  \\begin{axis}[ybar interval, ymax=" << T(double(p.max) * 1.1)  << ",ymin=0, minor y tick num = 3" << std::endl;
      //os << "    \\addplot coordinates { ";
    } else {
      os << "  \\begin{axis}[ymax=" << T(double(p.max) * 1.1)  << ",ymin=0" << std::endl;
    }
    if (p.xcaption.size()) {
      os << ", xlabel={" << p.xcaption << "}";
    }
    if (p.ycaption.size()) {
      os << ", ylabel={" << p.ycaption << "}";
    }

    os << "]" << std::endl;
    for (unsigned int pl = 0; pl < p.plots.size(); ++pl) {
      os << "    \\addplot";
      if (!p.histo) {
         os << "[" << "mark=" << p.plots_marks[pl];
         os <<", color=" << p.plots_colors[pl] << "]";
      }
      os << " coordinates { ";
      for (unsigned int i = 0; i < p.plots[pl]->size(); ++i) {
        os << "(" << i << "," << (*p.plots[pl])[i] << ") ";
      }
      if (p.histo) {
        // add one value because latex ignores the last one (i dont know why)
        os << "(" << p.plots[pl]->size() << ", " << p.max / 2 << ") ";
      }
      os << "};" << std::endl;
    }
    for (unsigned int l = 0; l < p.lines.size(); ++l) {
      os << "\\draw [" << p.lines[l].color  << (p.lines[l].dashed ? ", dashed" : "") 
         << "] ({rel axis cs:0,0}|-{axis cs:0," << p.lines[l].y 
         << "}) -- ({rel axis cs:1,0}|-{axis cs:" << p.plots[0]->size() << "," << p.lines[l].y 
         << "}) node [pos=0.5, below] {" << p.lines[l].caption << "};" << std::endl;
    }
    if (p.legend.size()) {
      os << "\\node[draw=black,thick,rounded corners=2pt,above right=2mm] at (0, 0) {%" << std::endl;
      os << "  \\begin{tabular}{@{}r@{ }l@{}}" << std::endl;
      for (unsigned int l = 0; l < p.legend.size(); ++l) {
        os << "    \\raisebox{2pt}{\\tikz{\\draw[" << p.legend[l].color  
           << "] (0,0) -- (5mm,0);}}&" << p.legend[l].text << "\\\\" << std::endl;
      }
      os << "  \\end{tabular}};" << std::endl;

    }
    for (unsigned int a = 0; a < p.arrows.size(); ++a) {
        os << "\\draw ["<< p.arrows[a].color <<   ", <->] ({axis cs:" << p.arrows[a].x + 0.5
           << "," << p.arrows[a].ymin << "}) -- ({axis cs:" << p.arrows[a].x + 0.5 
           << "," << p.arrows[a].ymax << "}) node [pos=0.5, left] { " << p.arrows[a].ymax - p.arrows[a].ymin
           << "};" << std::endl;
      
    }
    /*
    if ((double)variance > 0.01) {
        os << "\\node[red, draw, scale=1.2] at ({rel axis cs:0.5,0.25}) {variance : "
           << variance << "};" << std::endl;
    }
    */
    os << "  \\end{axis}" << std::endl;
    os << "\\end{tikzpicture}" << std::endl;
    os << "\\caption*{" << p.caption << "}" << std::endl;
    os << "\\end{minipage}" << std::endl;
    return os;
  }

  T max;
  std::string caption;
  bool histo;
  std::vector<const std::vector<T> *>plots;
  std::vector<std::string> plots_colors;
  std::vector<char> plots_marks;
  std::vector<HorizontalLine> lines;
  std::vector<Legend> legend;
  std::vector<Arrow> arrows;
  std::string xcaption;
  std::string ycaption;

};

class TexWriter {
public:
  TexWriter(const std::string &file): os(file.c_str()) { 
    if (!os.is_open()) {
      std::cerr << "Warning, cannot open : " << file << " for writing" << std::endl;
    }
  }
  ~TexWriter() {os.close();}
  
  template <typename T>
  void write_plot(Plot<T> &plot) { os << plot; }
private:
  std::ofstream os;
};

#endif

