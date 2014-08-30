
#include <map>

#include <TObject.h>


namespace rpwa {

	class dataMetadata : public TObject {
		friend class rootpwaDataFileWriter;

	  private:
		typedef std::pair<double, double> rangePairType;
		typedef std::map<std::string, rangePairType> binningMapType;

	  public:

		dataMetadata();
		~dataMetadata();

		const std::string& userString() const { return _userString; }
		const std::string& contentHash() const { return _contentHash; }
		const binningMapType& getBinningMap() const { return _binningMap; }
		const std::vector<std::string>& initalStateParticleNames() const { return _initialStateParticleNames; }
		const std::vector<std::string>& finalStateParticleNames() const { return _finalStateParticleNames; }

		std::ostream& print(std::ostream& out) const;


	  private:

		void setUserString(const std::string& userString) { _userString = userString; }
		void setContentHash(const std::string& contentHash) { _contentHash = contentHash; }
		void setInitialStateParticleNames(const std::vector<std::string>& initialStateParticleNames);
		void setFinalStateParticleNames(const std::vector<std::string>& finalStateParticleNames);

		void setBinningVariableLabels(const std::vector<std::string>& labels);
		void setBinningVariableRange(const std::string& label, const rangePairType& range);
		void setBinningMap(const binningMapType& binningMap);

		std::string _userString;
		std::string _contentHash;

		std::vector<std::string> _initialStateParticleNames;
		std::vector<std::string> _finalStateParticleNames;

		std::map<std::string, std::pair<double, double> > _binningMap;

		ClassDef(dataMetadata, 1);

	}; // class dataMetadata


	inline
	std::ostream&
	operator <<(std::ostream&          out,
	            const dataMetadata&    metadata)
	{
		return metadata.print(out);
	}

} // namespace rpwa
