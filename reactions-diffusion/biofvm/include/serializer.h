#pragma once
#include <memory>

namespace physicore::biofvm {

class microenvironment;

class serializer
{
public:
	virtual void serialize(const microenvironment& m) = 0;
};

using serializer_ptr = std::unique_ptr<serializer>;

} // namespace physicore::biofvm
