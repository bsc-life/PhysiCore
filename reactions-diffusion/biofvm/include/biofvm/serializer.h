#pragma once
#include <memory>

#include <biofvm/biofvm_export.h>
#include <common/types.h>

namespace physicore::reactions_diffusion::biofvm {

class microenvironment;

class BIOFVM_EXPORT serializer
{
public:
	virtual void serialize(const microenvironment& m, real_t current_time) = 0;

	virtual ~serializer() = default;
};

using serializer_ptr = std::unique_ptr<serializer>;

} // namespace physicore::reactions_diffusion::biofvm
