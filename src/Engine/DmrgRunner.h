#ifndef DMRGRUNNER_H
#define DMRGRUNNER_H
#include "BasisWithOperators.h"
#include "CmdLineOptions.hh"
#include "CrsMatrix.h"
#include "DmrgSolver.h"
#include "InputCheck.h"
#include "Introspect.hh"
#include "LeftRightSuper.h"
#include "MatrixVectorKron/MatrixVectorKron.h"
#include "MatrixVectorOnTheFly.h"
#include "MatrixVectorStored.h"
#include "ModelHelperLocal.h"
#include "ModelSelector.h"
#include "ParametersDmrgSolver.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"
#include "Qn.h"
#include "RedirectOutput.hh"
#include "RunFinished.h"
#include "SuperGeometry.h"
#include "VectorWithOffset.h"

namespace Dmrg {

template <typename RealType> class DmrgRunner {

public:

	using InputNgType              = PsimagLite::InputNg<InputCheck>;
	using ParametersDmrgSolverType = ParametersDmrgSolver<RealType, InputNgType::Readable, Qn>;
	using ApplicationType          = PsimagLite::PsiApp;

	DmrgRunner(const ApplicationType& application)
	    : application_(application)
	{ }

	void doOneRun(const PsimagLite::String& data, const CmdLineOptions& cmd_line) const;

	void doOneRun(const PsimagLite::String&   data,
	              const CmdLineOptions&       cmd_line,
	              const OptionsForIntrospect& op_options) const;

	const ApplicationType& application() const { return application_; }

private:

	template <typename ComplexOrRealType>
	void doOneRun2(const ParametersDmrgSolverType& dmrgSolverParams,
	               const OptionsForIntrospect&     op_options,
	               InputNgType::Readable&          io) const;

	template <typename MatrixVectorType>
	void doOneRun3(const ParametersDmrgSolverType& dmrgSolverParams,
	               const OptionsForIntrospect&     op_options,
	               InputNgType::Readable&          io) const;

	template <typename MatrixVectorType, typename VectorWithOffsetType>
	void doOneRun4(const ParametersDmrgSolverType& dmrgSolverParams,
	               const OptionsForIntrospect&     op_options,
	               InputNgType::Readable&          io) const;

	void dealWithConsoleOutput(const ParametersDmrgSolverType& dmrgSolverParams,
	                           const std::string&              logfile,
	                           bool                            unbuffered) const;

	void adjustDmrgSolverParams(ParametersDmrgSolverType& dmrgSolverParams,
	                            const std::string&        insitu,
	                            SizeType                  precision,
	                            SizeType                  threadsInCmd) const;

	static void setCorrectThreading(const ParametersDmrgSolverType& dmrgSolverParams,
	                                InputNgType::Readable&          io);

	static bool endDueToClobber(bool is_no_clobber_set, const std::string& logfile);

	const ApplicationType& application_;
};
}
#endif // DMRGRUNNER_H
