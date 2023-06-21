#include "cellmlmathml.h"
#include <limits>
#include <symengine/eval_double.h>
#include <symengine/parser.h>
#include <symengine/printers.h>



std::string replaceAll(std::string str, const std::string &from,
                       const std::string &to) {
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos +=
        to.length(); // Handles case where 'to' is a substring of 'from'
  }
  return str;
}

namespace SymEngine {

std::vector<std::string> init_mathml_printer_names() {
  std::vector<std::string> names = init_str_printer_names();
  names[SYMENGINE_ASIN] = "arcsin";
  names[SYMENGINE_ACOS] = "arccos";
  names[SYMENGINE_ASEC] = "arcsec";
  names[SYMENGINE_ACSC] = "arccsc";
  names[SYMENGINE_ATAN] = "arctan";
  names[SYMENGINE_ACOT] = "arccot";
  names[SYMENGINE_ASINH] = "arcsinh";
  names[SYMENGINE_ACSCH] = "arccsch";
  names[SYMENGINE_ACOSH] = "arccosh";
  names[SYMENGINE_ATANH] = "arctanh";
  names[SYMENGINE_ACOTH] = "arccoth";
  names[SYMENGINE_ASECH] = "arcsech";
  return names;
}

void CellMLMathMLPrinter::bvisit(const Basic &x) {
  throw SymEngineException("Error: not supported");
}

void CellMLMathMLPrinter::bvisit(const Symbol &x) {
  s << "<ci>" << x.get_name() << "</ci>";
}

void CellMLMathMLPrinter::bvisit(const Integer &x) {
  s << "<cn cellml:units=\"dimensionless\">" << x.as_integer_class() << "</cn>";
}

void CellMLMathMLPrinter::bvisit(const Rational &x) {
  const auto &rational = x.as_rational_class();
  s << "<cn cellml:units=\"dimensionless\">" << get_num(rational) << "<sep/>"
    << get_den(rational) << "</cn>";
}

void CellMLMathMLPrinter::bvisit(const RealDouble &x) {
  s << "<cn cellml:units=\"dimensionless\">" << x << "</cn>";
}

#ifdef HAVE_SYMENGINE_MPFR
void CellMLMathMLPrinter::bvisit(const RealMPFR &x) {
  // TODO: Use bigfloat here
  s << "<cn cellml:units=\"dimensionless\">" << x << "</cn>";
}
#endif

void CellMLMathMLPrinter::bvisit(const ComplexBase &x) {
  s << "<apply><csymbol cd=\"nums1\">complex_cartesian</csymbol>";
  x.real_part()->accept(*this);
  x.imaginary_part()->accept(*this);
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const Interval &x) {
  s << "<interval closure=";
  if (x.get_left_open()) {
    if (x.get_right_open()) {
      s << "\"open\">";
    } else {
      s << "\"open-closed\">";
    }
  } else {
    if (x.get_right_open()) {
      s << "\"closed-open\">";
    } else {
      s << "\"closed\">";
    }
  }
  x.get_start()->accept(*this);
  x.get_end()->accept(*this);
  s << "</interval>";
}

void CellMLMathMLPrinter::bvisit(const Piecewise &x) {
  s << "<piecewise>";
  const auto &equations = x.get_vec();
  for (const auto &equation : equations) {
    s << "<piece>";
    equation.first->accept(*this);
    equation.second->accept(*this);
    s << "</piece>";
  }
  s << "</piecewise>";
}

void CellMLMathMLPrinter::bvisit(const EmptySet &x) { s << "<emptyset/>"; }

void CellMLMathMLPrinter::bvisit(const Complexes &x) { s << "<complexes/>"; }

void CellMLMathMLPrinter::bvisit(const Reals &x) { s << "<reals/>"; }

void CellMLMathMLPrinter::bvisit(const Rationals &x) { s << "<rationals/>"; }

void CellMLMathMLPrinter::bvisit(const Integers &x) { s << "<integers/>"; }

void CellMLMathMLPrinter::bvisit(const FiniteSet &x) {
  s << "<set>";
  const auto &args = x.get_args();
  for (const auto &arg : args) {
    arg->accept(*this);
  }
  s << "</set>";
}

void CellMLMathMLPrinter::bvisit(const ConditionSet &x) {
  s << "<set><bvar>";
  x.get_symbol()->accept(*this);
  s << "</bvar><condition>";
  x.get_condition()->accept(*this);
  s << "</condition>";
  x.get_symbol()->accept(*this);
  s << "</set>";
}

void CellMLMathMLPrinter::bvisit(const Contains &x) {
  s << "<apply><in/>";
  x.get_expr()->accept(*this);
  x.get_set()->accept(*this);
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const BooleanAtom &x) {
  if (x.get_val()) {
    s << "<true/>";
  } else {
    s << "<false/>";
  }
}

void CellMLMathMLPrinter::bvisit(const And &x) {
  s << "<apply><and/>";
  const auto &conditions = x.get_args();
  for (const auto &condition : conditions) {
    condition->accept(*this);
  }
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const Or &x) {
  s << "<apply><or/>";
  const auto &conditions = x.get_args();
  for (const auto &condition : conditions) {
    condition->accept(*this);
  }
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const Xor &x) {
  s << "<apply><xor/>";
  const auto &conditions = x.get_args();
  for (const auto &condition : conditions) {
    condition->accept(*this);
  }
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const Not &x) {
  s << "<apply><not/>";
  x.get_arg()->accept(*this);
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const Union &x) {
  s << "<apply><union/>";
  const auto &sets = x.get_args();
  for (const auto &set : sets) {
    set->accept(*this);
  }
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const Complement &x) {
  s << "<apply><setdiff/>";
  x.get_universe()->accept(*this);
  x.get_container()->accept(*this);
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const ImageSet &x) {
  s << "<set><bvar>";
  x.get_expr()->accept(*this);
  s << "</bvar><condition><apply><in/>";
  x.get_symbol()->accept(*this);
  x.get_baseset()->accept(*this);
  s << "</apply></condition>";
  x.get_symbol()->accept(*this);
  s << "</set>";
}

void CellMLMathMLPrinter::bvisit(const Add &x) {
  s << "<apply><plus/>";
  auto args = x.get_args();
  for (auto arg : args) {
    arg->accept(*this);
  }
  s << "</apply>";
}

std::map<std::string,std::string> CellMLMathMLPrinter::mulexpressionMap = {};

void CellMLMathMLPrinter::bvisit(const Mul &x) {
  //Reuse if encountered apriori
  std::string x_expr = x.__str__();
  if(mulexpressionMap.find(x_expr)!=mulexpressionMap.end()){
    //std::cout<<"Resuing "<<mulexpressionMap[x_expr]<< " for "<<x_expr<<std::endl;
    s << mulexpressionMap[x_expr];
    return;
  }

  // Store current state as s gets altered by children
  std::string current = s.str();
  s.str("");
  s.clear();
  std::vector<SymEngine::TypeID> atypes;
  std::vector<std::string> mathml;
  std::vector<bool> ignore;
  std::vector<bool> denom;
  bool numAvailable = false;
  unsigned int divisionOp = 0;
  unsigned int leftCount = 0;
  auto args = x.get_args();
  size_t negOne = -1;
  for (int i = 0; i < args.size(); i++) {
    const auto &arg = args[i];
    auto code = arg->get_type_code();
    // Handle negative symbols
    if (code == SymEngine::SYMENGINE_INTEGER &&
        negativeOne->compare(*arg) == 0) {
      negOne = i;
      numAvailable = true;
    }
    if (code == SymEngine::SYMENGINE_SYMBOL) {
      numAvailable = true;
    }
    atypes.push_back(code);
    arg->accept(*this);
    mathml.push_back(s.str());
    ignore.push_back(false);
    denom.push_back(false);

    s.str("");
    s.clear();
  }
  // Check if there is any negative 1
  if (negOne != -1) {
    // Find first symbol and make it -symbol
    for (int a = 0; a < atypes.size(); a++) {
      if (atypes[a] == SymEngine::SYMENGINE_SYMBOL) {
        ignore[negOne] = true;
        mathml[a] = "<apply><minus/>" + mathml[a] + "</apply>";
        break;
      }
    }
  }
  // Handle division
  for (int a = 0; a < atypes.size(); a++) {
    if (atypes[a] == SymEngine::SYMENGINE_POW) {
      if (mathml[a].find(
              "<apply><divide/><cn cellml:units=\"dimensionless\">1</cn>") !=
          std::string::npos) {
        std::string ax = replaceAll(
            mathml[a],
            "<apply><divide/><cn cellml:units=\"dimensionless\">1</cn>", "");
        ax = replaceAll(ax, "</apply>", "");
        mathml[a] = ax;
        denom[a] = true;
        divisionOp++;
      }else if((mathml[a].find("<apply><power/>") != std::string::npos) && 
        (mathml[a].find(
              "<cn cellml:units=\"dimensionless\">-1</cn></apply>") !=
          std::string::npos)){
        std::string ax = replaceAll(mathml[a],"<apply><power/>", "");
        ax = replaceAll(ax,"<cn cellml:units=\"dimensionless\">-1</cn></apply>", "");
        mathml[a] = ax;
        denom[a] = true;
        divisionOp++;
      }
    }
  }
  s.str("");
  s.clear();
  // Handle divisions
  if (divisionOp > 0) {
    std::ostringstream ss;
    if (divisionOp > 1) {
      ss << "<apply><times/>";
      for (int a = 0; a < atypes.size(); a++) {
        if (denom[a]) {
          ss << mathml[a];
          ignore[a] = true;
        }
      }
      ss << "</apply>";
    } else {
      for (int a = 0; a < atypes.size(); a++) {
        if (denom[a]) {
          ss << mathml[a];
          ignore[a] = true;
        }
      }
    }
    // Find a numerator
    if (numAvailable) {
      for (int a = atypes.size() - 1; a > -1; a--) {
        if (!denom[a] && !ignore[a]) {
          s << "<apply><divide/>" << mathml[a];
          ignore[a] = true;
          break;
        }
      }
    } else {
      // If no numerators are available use 1
      s << "<apply><divide/><cn cellml:units=\"dimensionless\">1.0</cn>";
    }
    // Check if any terms appearing in the numerator have not been included
    leftCount = 0;
    for (int a = 0; a < atypes.size(); a++) {
      if (!ignore[a]) {
        leftCount++;
      }
    }
    if (leftCount > 0) {
      s << "<apply><times/>";
      for (int a = 0; a < atypes.size(); a++) {
        if (!ignore[a]) {
          s << mathml[a];
        }
      }
    }
    s << ss.str();
    // Close apply for divide
    s << "</apply>";
  } else {
    // Check if more than one term is left to use times
    leftCount = 0;
    for (int a = 0; a < atypes.size(); a++) {
      if (!ignore[a]) {
        leftCount++;
      }
    }
    // If only one term is left - typically as -x (as we are in Mul), it's
    // mathml is <apply><minus/>x</apply>
    if (leftCount > 1) {
      s << "<apply><times/>";
    }
    for (int a = 0; a < atypes.size(); a++) {
      if (!ignore[a]) {
        s << mathml[a];
      }
    }
    if (leftCount > 1) {
      s << "</apply>";
    }
  }
  // s has xml for this op
  std::string mulstring = s.str();
  mulexpressionMap[x_expr] = mulstring;
  s.str("");
  s.clear();

  
  // Add the previous state to s and return
  s << current << mulstring;

  // Default code from symengine - left for reference
  //   auto args = x.get_args();
  //   for (auto arg : args) {
  //     arg->accept(*this);
  //   }
  //   s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const Pow &x) {
  // Compare is similar to string compare
  const auto &exp = x.get_exp();
  if (exp->get_type_code() == SymEngine::SYMENGINE_CONSTANT &&
      exp->compare(*negativeOne) == 0) {
    s << "<apply><divide/><cn cellml:units=\"dimensionless\">1</cn>";
    x.get_base()->accept(*this);
    s << "</apply>";
  } else {
    s << "<apply><power/>";
    x.get_base()->accept(*this);
    x.get_exp()->accept(*this);
    s << "</apply>";
  }
}

void CellMLMathMLPrinter::bvisit(const Constant &x) {
  s << "<";
  if (eq(x, *pi)) {
    s << "pi/";
  } else if (eq(x, *E)) {
    s << "exponentiale/";
  } else if (eq(x, *EulerGamma)) {
    s << "eulergamma/";
  } else {
    s << "cn cellml:units=\"dimensionless\">" << eval_double(x) << "</cn";
  }
  s << ">";
}

void CellMLMathMLPrinter::bvisit(const Function &x) {
  static const std::vector<std::string> names_ = init_mathml_printer_names();
  s << "<apply>";
  s << "<" << names_[x.get_type_code()] << "/>";
  const auto &args = x.get_args();
  for (const auto &arg : args) {
    arg->accept(*this);
  }
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const UnevaluatedExpr &x) {
  apply(*x.get_arg());
}

void CellMLMathMLPrinter::bvisit(const FunctionSymbol &x) {
  s << "<apply><ci>" << x.get_name() << "</ci>";
  const auto &args = x.get_args();
  for (const auto &arg : args) {
    arg->accept(*this);
  }
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const Equality &x) {
  s << "<apply><eq/>";
  x.get_arg1()->accept(*this);
  x.get_arg2()->accept(*this);
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const Unequality &x) {
  s << "<apply><neq/>";
  x.get_arg1()->accept(*this);
  x.get_arg2()->accept(*this);
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const LessThan &x) {
  s << "<apply><leq/>";
  x.get_arg1()->accept(*this);
  x.get_arg2()->accept(*this);
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const StrictLessThan &x) {
  s << "<apply><lt/>";
  x.get_arg1()->accept(*this);
  x.get_arg2()->accept(*this);
  s << "</apply>";
}

void CellMLMathMLPrinter::bvisit(const Derivative &x) {
  s << "<apply><partialdiff/><bvar>";
  for (const auto &elem : x.get_symbols()) {
    elem->accept(*this);
  }
  s << "</bvar>";
  x.get_arg()->accept(*this);
  s << "</apply>";
}

std::string CellMLMathMLPrinter::apply(const Basic &b) {
  b.accept(*this);
  return s.str();
}

std::string cellmlmathml(const Basic &x) {
  CellMLMathMLPrinter m;
  return m.apply(x);
}
} // namespace SymEngine
