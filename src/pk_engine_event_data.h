// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains event-data ownership, views, ADDL guards, and structural signatures.

struct OwnedEventTable {
  std::unordered_map<std::string, std::vector<double>> columns;
  int rows = 0;

  OwnedEventTable(const Rcpp::DataFrame& source, int start, int length)
      : rows(length) {
    if (start < 0 || length < 0 || start > source.nrows() ||
        length > source.nrows() - start) {
      throw std::out_of_range("Native event-data copy is outside its source table.");
    }
    Rcpp::CharacterVector names = source.names();
    columns.reserve(static_cast<std::size_t>(names.size()));
    for (R_xlen_t column_index = 0; column_index < names.size(); ++column_index) {
      const std::string name = Rcpp::as<std::string>(names[column_index]);
      Rcpp::RObject column = source[column_index];
      const SEXPTYPE type = static_cast<SEXPTYPE>(TYPEOF(column));
      if (type != REALSXP && type != INTSXP && type != LGLSXP) continue;
      std::vector<double> values(static_cast<std::size_t>(length));
      if (type == REALSXP) {
        const double* input = REAL_RO(column);
        std::copy(input + start, input + start + length, values.begin());
      } else {
        const int* input = INTEGER_RO(column);
        for (int row = 0; row < length; ++row) {
          const int value = input[start + row];
          values[static_cast<std::size_t>(row)] =
            value == NA_INTEGER ? NA_REAL : static_cast<double>(value);
        }
      }
      columns.emplace(name, std::move(values));
    }
  }

  bool contains(const std::string& name) const {
    return columns.find(name) != columns.end();
  }
};

// A non-owning row-range view over one canonical R event data frame.  R data
// frame subsetting duplicates every selected column; this view instead keeps
// the parent columns alive and translates local subject rows to their original
// offsets.  Full-population execution is represented by the same type with a
// zero offset, so numerical code has one data-access path.  The alternative
// OwnedEventTable constructor is used by guarded native value-only workers.
class EventDataView {
 public:
  class Column {
   public:
    class Iterator {
     public:
      Iterator(const Column* column, int row) : column_(column), row_(row) {}
      double operator*() const { return (*column_)[row_]; }
      Iterator& operator++() { ++row_; return *this; }
      bool operator!=(const Iterator& other) const { return row_ != other.row_; }
     private:
      const Column* column_;
      int row_;
    };

    Column(const EventDataView* view, std::string name,
           double fallback = std::numeric_limits<double>::quiet_NaN(),
           bool allow_missing = false);
    double operator[](int row) const;
    int size() const;
    Iterator begin() const { return Iterator(this, 0); }
    Iterator end() const { return Iterator(this, size()); }

   private:
    std::string name_;
    SEXP column_ = R_NilValue;
    const double* real_data_ = nullptr;
    const int* integer_data_ = nullptr;
    const double* native_data_ = nullptr;
    std::shared_ptr<const OwnedEventTable> native_owner_;
    SEXPTYPE type_ = NILSXP;
    int start_ = 0;
    int rows_ = 0;
    bool subject_view_ = false;
    bool missing_ = false;
    double fallback_;
    bool allow_missing_;
  };

  EventDataView(const Rcpp::DataFrame& data)
      : data_(std::make_shared<Rcpp::DataFrame>(data)), start_(0),
        rows_(data.nrows()), subject_view_(false) {}

  EventDataView(const Rcpp::DataFrame& data, int start, int rows,
                bool subject_view = true)
      : data_(std::make_shared<Rcpp::DataFrame>(data)), start_(start),
        rows_(rows), subject_view_(subject_view) {
    if (start_ < 0 || rows_ < 0 || start_ > data.nrows() ||
        rows_ > data.nrows() - start_) {
      throw std::out_of_range("Event-data row view is outside its parent table.");
    }
  }

  explicit EventDataView(std::shared_ptr<const OwnedEventTable> data,
                         bool subject_view = true)
      // Do not call DataFrame::create() here: native views are constructed in
      // optimized subject workers, where any R allocation/PROTECT operation is
      // forbidden.  The default wrapper remains NIL and is never accessed
      // while native_ is present.
      : data_(nullptr), native_(std::move(data)), start_(0),
        rows_(native_ ? native_->rows : 0), subject_view_(subject_view) {
    if (!native_ || rows_ < 1) {
      throw std::invalid_argument("Native event-data view is empty.");
    }
  }

  int nrows() const { return rows_; }
  int start() const { return start_; }
  bool is_subject_view() const { return subject_view_; }
  bool containsElementNamed(const char* name) const {
    return native_ ? native_->contains(name) :
      data_->containsElementNamed(name);
  }
  const Rcpp::DataFrame& parent() const {
    if (native_) {
      throw std::logic_error("A native event-data view has no R parent.");
    }
    return *data_;
  }

  int source_row(int row) const {
    if (row < 0 || row >= rows_) {
      throw std::out_of_range("Event-data row is outside the active view.");
    }
    return start_ + row;
  }

  Rcpp::RObject column(const std::string& name) const {
    if (native_) {
      throw std::logic_error("Native event-data columns are not R objects.");
    }
    if (!containsElementNamed(name.c_str())) {
      throw std::invalid_argument(
        "PRED input '" + name + "' is not present in the dataset.");
    }
    auto cached = columns_.find(name);
    if (cached != columns_.end()) return cached->second;
    Rcpp::RObject value = (*data_)[name];
    columns_.emplace(name, value);
    return value;
  }
  Column values(const std::string& name) const { return Column(this, name); }
  Column values(const std::string& name, double fallback) const {
    return Column(this, name, fallback, true);
  }

 private:
  std::shared_ptr<Rcpp::DataFrame> data_;
  std::shared_ptr<const OwnedEventTable> native_;
  int start_ = 0;
  int rows_ = 0;
  bool subject_view_ = false;
  mutable std::unordered_map<std::string, Rcpp::RObject> columns_;
};

inline EventDataView::Column::Column(
    const EventDataView* view, std::string name, double fallback,
    bool allow_missing)
    : name_(std::move(name)), start_(view->start()),
      rows_(view->nrows()), subject_view_(view->is_subject_view()),
      missing_(!view->containsElementNamed(name_.c_str())),
      fallback_(fallback), allow_missing_(allow_missing) {
  if (view->native_) {
    native_owner_ = view->native_;
    auto found = native_owner_->columns.find(name_);
    missing_ = found == native_owner_->columns.end();
    if (missing_ && !allow_missing_) {
      throw std::invalid_argument(
        "PRED input '" + name_ + "' is not present in the native dataset.");
    }
    if (!missing_) native_data_ = found->second.data();
    return;
  }
  if (!missing_) column_ = view->column(name_);
  if (missing_ && !allow_missing_) {
    throw std::invalid_argument(
      "PRED input '" + name_ + "' is not present in the dataset.");
  }
  if (!missing_) {
    type_ = static_cast<SEXPTYPE>(TYPEOF(column_));
    // Resolve ALTREP and cache a read-only address once, on R's main thread.
    // The retained RObject protects the backing vector for this view's whole
    // lifetime. Subsequent row access therefore performs no Rcpp wrapper
    // construction and no R API call.
    if (type_ == REALSXP) {
      real_data_ = REAL_RO(column_);
    } else if (type_ == INTSXP || type_ == LGLSXP) {
      integer_data_ = INTEGER_RO(column_);
    }
  }
}

inline double EventDataView::Column::operator[](int row) const {
  if (row < 0 || row >= rows_) {
    throw std::out_of_range("Event-data row is outside the active column view.");
  }
  if (subject_view_ && name_ == ".ID_INDEX") {
    return 1.0;
  }
  if (missing_) return fallback_;
  const int source = start_ + row;
  if (native_data_) return native_data_[source];
  if (type_ == REALSXP) return real_data_[source];
  if (type_ == INTSXP || type_ == LGLSXP) {
    const int value = integer_data_[source];
    return value == NA_INTEGER ? NA_REAL : static_cast<double>(value);
  }
  throw std::invalid_argument("Event-data column '" + name_ + "' must be numeric.");
}

inline int EventDataView::Column::size() const { return rows_; }

// Retains the canonical data frame once and describes every contiguous subject
// with two integers.  The R-side environment owns this external pointer, while
// small nm_subject_view descriptors select a subject without copying columns.
struct NativeSubjectStore {
  Rcpp::DataFrame data;
  std::vector<int> starts;
  std::vector<int> lengths;

  NativeSubjectStore(const Rcpp::DataFrame& source,
                     const Rcpp::IntegerVector& starts_input,
                     const Rcpp::IntegerVector& lengths_input)
      : data(source), starts(Rcpp::as<std::vector<int>>(starts_input)),
        lengths(Rcpp::as<std::vector<int>>(lengths_input)) {
    if (starts.empty() || starts.size() != lengths.size()) {
      throw std::invalid_argument("Native subject ranges are invalid.");
    }
    for (std::size_t index = 0; index < starts.size(); ++index) {
      --starts[index];
      if (starts[index] < 0 || lengths[index] < 1 ||
          starts[index] > data.nrows() ||
          lengths[index] > data.nrows() - starts[index]) {
        throw std::out_of_range("A native subject range is outside the event table.");
      }
      if (index && starts[index] != starts[index - 1] + lengths[index - 1]) {
        throw std::invalid_argument("Native subject ranges must be contiguous.");
      }
    }
  }

  EventDataView view(int subject) const {
    if (subject < 0 || subject >= static_cast<int>(starts.size())) {
      throw std::out_of_range("Subject view is outside the native data store.");
    }
    return EventDataView(data, starts[static_cast<std::size_t>(subject)],
                         lengths[static_cast<std::size_t>(subject)]);
  }
};

double data_value(const EventDataView& data, const std::string& name, int row);

// A compiled set of MU equations.  MU expressions use the same LibeRtAD IR
// as $PK/$PRED, which keeps model-language semantics in one place while
// allowing subject-level covariates to be read directly from NativeSubjectStore
// row views.  Each program output is mapped to its (zero-based) ETA column.
struct NativeMuProgram {
  std::shared_ptr<const libertad::Program> program;
  std::vector<CompiledInputBinding> inputs;
  std::vector<std::size_t> outputs;
  std::vector<int> eta;

  NativeMuProgram(const Rcpp::List& ir, const Rcpp::IntegerVector& eta_input)
      : program(std::make_shared<const libertad::Program>(ir)),
        inputs(compile_input_bindings(program->input_names)),
        eta(Rcpp::as<std::vector<int>>(eta_input)) {
    outputs.resize(program->output_names.size());
    std::iota(outputs.begin(), outputs.end(), 0U);
    if (outputs.size() != eta.size() || outputs.empty()) {
      throw std::invalid_argument(
        "Compiled MU outputs and ETA destinations must have equal non-zero length.");
    }
    for (int& index : eta) {
      --index;
      if (index < 0) {
        throw std::invalid_argument("A compiled MU ETA destination is invalid.");
      }
    }
    for (const CompiledInputBinding& binding : inputs) {
      if (binding.eta >= 0 || binding.sigma >= 0 || binding.state >= 0 ||
          binding.zero_error) {
        throw std::invalid_argument(
          "MU expressions may depend only on THETAs and subject-level data.");
      }
    }
  }
};

inline std::vector<double> evaluate_mu_program(
    const NativeMuProgram& compiled, const EventDataView& data,
    const Rcpp::NumericVector& theta) {
  if (!data.nrows()) {
    throw std::invalid_argument("A MU program cannot evaluate an empty subject view.");
  }
  std::vector<double> inputs(compiled.inputs.size(), 0.0);
  for (std::size_t index = 0; index < compiled.inputs.size(); ++index) {
    const CompiledInputBinding& binding = compiled.inputs[index];
    if (binding.theta >= 0) {
      if (binding.theta >= theta.size()) {
        throw std::out_of_range("MU THETA index exceeds supplied values.");
      }
      inputs[index] = theta[binding.theta];
    } else {
      inputs[index] = data_value(data, binding.name, 0);
    }
    if (!std::isfinite(inputs[index])) {
      throw std::domain_error(
        "MU input '" + binding.name + "' is non-finite for a subject.");
    }
  }
  const std::vector<double> result =
    compiled.program->eval_outputs(inputs, compiled.outputs);
  for (double value : result) {
    if (!std::isfinite(value)) {
      throw std::domain_error("A compiled MU equation returned a non-finite value.");
    }
  }
  return result;
}

EventDataView event_data_view(SEXP input) {
  if (Rf_inherits(input, "data.frame")) {
    return EventDataView(Rcpp::as<Rcpp::DataFrame>(input));
  }
  if (!Rf_inherits(input, "nm_subject_view") || TYPEOF(input) != VECSXP) {
    throw std::invalid_argument(
      "Native execution requires a data frame or nm_subject_view.");
  }
  Rcpp::List descriptor(input);
  if (!descriptor.containsElementNamed("pointer") ||
      !descriptor.containsElementNamed("subject")) {
    throw std::invalid_argument("An nm_subject_view descriptor is incomplete.");
  }
  Rcpp::XPtr<NativeSubjectStore> store(descriptor["pointer"]);
  return store->view(Rcpp::as<int>(descriptor["subject"]) - 1);
}

double data_value(const EventDataView& data, const std::string& name, int row) {
  if (data.is_subject_view() && name == ".ID_INDEX") {
    data.source_row(row);
    return 1.0;
  }
  if (!data.containsElementNamed(name.c_str())) {
    throw std::invalid_argument("PRED input '" + name + "' is not present in the dataset.");
  }
  return data.values(name)[row];
}

void require_materialized_addl(const EventDataView& data) {
  if (!data.containsElementNamed("ADDL")) return;
  for (int row = 0; row < data.nrows(); ++row) {
    const double addl = data_value(data, "ADDL", row);
    if (!std::isfinite(addl) || addl != 0.0) {
      throw std::invalid_argument(
        "Native execution requires ADDL/II doses to be materialized by "
        "nm_dataset(); a non-zero or invalid ADDL value reached C++.");
    }
  }
}

inline void fnv_mix(std::uint64_t& hash, const void* source, std::size_t size,
                    std::uint64_t prime) {
  const unsigned char* bytes = static_cast<const unsigned char*>(source);
  for (std::size_t index = 0; index < size; ++index) {
    hash ^= static_cast<std::uint64_t>(bytes[index]);
    hash *= prime;
  }
}

inline void fnv_mix_string(std::uint64_t& first, std::uint64_t& second,
                           const std::string& value) {
  fnv_mix(first, value.data(), value.size(), UINT64_C(1099511628211));
  fnv_mix(second, value.data(), value.size(), UINT64_C(14029467366897019727));
  const unsigned char boundary = 0xffU;
  fnv_mix(first, &boundary, 1U, UINT64_C(1099511628211));
  fnv_mix(second, &boundary, 1U, UINT64_C(14029467366897019727));
}

std::string event_data_signature(
    const EventDataView& data,
    const std::unordered_map<std::string, bool>& ignored,
    bool include_fo_layout) {
  std::uint64_t first = UINT64_C(1469598103934665603);
  std::uint64_t second = UINT64_C(7809847782465536322);
  Rcpp::CharacterVector names = data.parent().names();
  for (R_xlen_t column_index = 0; column_index < names.size(); ++column_index) {
    const std::string name = Rcpp::as<std::string>(names[column_index]);
    if (ignored.find(name) != ignored.end()) continue;
    fnv_mix_string(first, second, name);
    Rcpp::RObject column = data.parent()[column_index];
    const int type = TYPEOF(column);
    fnv_mix(first, &type, sizeof(type), UINT64_C(1099511628211));
    fnv_mix(second, &type, sizeof(type), UINT64_C(14029467366897019727));
    for (int row = 0; row < data.nrows(); ++row) {
      const int source = data.source_row(row);
      if (data.is_subject_view() && name == ".STRUCT_ID_INDEX") {
        const int current = Rcpp::IntegerVector(column)[source];
        const unsigned char boundary = row == 0 || current !=
          Rcpp::IntegerVector(column)[data.source_row(row - 1)];
        fnv_mix(first, &boundary, sizeof(boundary), UINT64_C(1099511628211));
        fnv_mix(second, &boundary, sizeof(boundary), UINT64_C(14029467366897019727));
        continue;
      }
      if (type == REALSXP) {
        const double value = Rcpp::NumericVector(column)[source];
        fnv_mix(first, &value, sizeof(value), UINT64_C(1099511628211));
        fnv_mix(second, &value, sizeof(value), UINT64_C(14029467366897019727));
      } else if (type == INTSXP || type == LGLSXP) {
        const int value = Rcpp::IntegerVector(column)[source];
        fnv_mix(first, &value, sizeof(value), UINT64_C(1099511628211));
        fnv_mix(second, &value, sizeof(value), UINT64_C(14029467366897019727));
      } else if (type == STRSXP) {
        SEXP value = STRING_ELT(column, source);
        fnv_mix_string(first, second,
          value == NA_STRING ? std::string("<NA>") : std::string(CHAR(value)));
      } else {
        throw std::invalid_argument(
          "Subject structural signatures require atomic event-data columns.");
      }
    }
  }
  if (include_fo_layout) {
    fnv_mix_string(first, second, "<FO_LAYOUT>");
    const bool has_dvid = data.containsElementNamed("DVID");
    for (int row = 0; row < data.nrows(); ++row) {
      const unsigned char observed =
        data_value(data, "EVID", row) == 0.0 &&
        data_value(data, "MDV", row) == 0.0 &&
        std::isfinite(data_value(data, "DV", row));
      const int dvid = has_dvid ?
        std::max(1, static_cast<int>(data_value(data, "DVID", row))) : 1;
      fnv_mix(first, &observed, sizeof(observed), UINT64_C(1099511628211));
      fnv_mix(second, &observed, sizeof(observed), UINT64_C(14029467366897019727));
      fnv_mix(first, &dvid, sizeof(dvid), UINT64_C(1099511628211));
      fnv_mix(second, &dvid, sizeof(dvid), UINT64_C(14029467366897019727));
    }
  }
  std::ostringstream output;
  output << std::hex << std::setfill('0') << std::setw(16) << first
         << std::setw(16) << second;
  return output.str();
}

