#pragma once
#include <memory>

#include <common/types.h>

namespace physicore::biofvm {

class microenvironment;

class serializer
{
public:
	virtual void serialize(const microenvironment& m, real_t current_time) = 0;

	virtual ~serializer() = default;
};

using serializer_ptr = std::unique_ptr<serializer>;

} // namespace physicore::biofvm
